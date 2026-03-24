import os
import subprocess

# Geometry creation follows two paths:
#   1) use an explicit .geo/.step file if one exists
#   2) otherwise, parse a limited subset of MCNP surface cards and build a
#      simple equivalent geometry through the Gmsh Python API


def _normalize_text_file_to_lf(path):
    # OpenSn expects Unix newlines in text mesh files.
    # Normalize CRLF/CR to LF after Gmsh writes the file.
    with open(path, "rb") as f:
        raw = f.read()

    # Skip binary files.
    if b"\x00" in raw:
        return

    normalized = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    if normalized == raw:
        return

    with open(path, "wb") as f:
        f.write(normalized)


class Plane:
    def __init__(self, axis, value):
        self.axis = axis
        self.value = value


class CylinderSurf:
    def __init__(self, axis, x0, y0, z0, dx, dy, dz, R):
        self.axis = axis
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.R = R


class SphereSurf:
    def __init__(self, x0, y0, z0, R):
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.R = R


def _detect_geometry_kind_from_name(path):
    stem = os.path.splitext(os.path.basename(path))[0].upper()
    for suffix in ("-IN", "-SL", "-CY", "-SP", "-ISLC"):
        if stem.endswith(suffix):
            return suffix[1:]
    return "AUTO"


def _parse_mcnp_surfaces(mcnp_input):
    # Extract just enough structured geometry information from the MCNP input
    # to build simple slab/sphere/cylinder meshes.
    px_list = []
    py_list = []
    pz_list = []
    cyl_list = []
    sph_list = []

    with open(mcnp_input, "r") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if "$" in line:
                line = line.split("$", 1)[0]
            if not line.strip():
                continue
            if line[0].lower() == "c":
                continue

            tokens = line.split()
            if len(tokens) < 3:
                continue
            if not tokens[0].lstrip("+-").isdigit():
                continue

            # Accept MCNP variants like "c/z" as well as "cz".
            mnemonic = tokens[1].lower().replace("/", "")

            # Planes PX/PY/PZ
            if mnemonic in ("px", "py", "pz"):
                try:
                    val = float(tokens[2])
                except ValueError:
                    continue
                axis = mnemonic[1]
                plane = Plane(axis=axis, value=val)
                if axis == "x":
                    px_list.append(plane)
                elif axis == "y":
                    py_list.append(plane)
                else:
                    pz_list.append(plane)
                continue

            # Sphere at origin: SO R
            if mnemonic == "so":
                try:
                    R = float(tokens[2])
                except ValueError:
                    continue
                if R > 0.0:
                    sph_list.append(SphereSurf(0.0, 0.0, 0.0, R))
                continue

            # General sphere: S/SPH x0 y0 z0 R
            if mnemonic in ("s", "sph") and len(tokens) >= 6:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    z0 = float(tokens[4])
                    R = float(tokens[5])
                except ValueError:
                    continue
                if R > 0.0:
                    sph_list.append(SphereSurf(x0, y0, z0, R))
                continue

            # Infinite cylinders: CZ/CX/CY
            # Accept both short forms (e.g. "cz R") and shifted forms
            # (e.g. "cz x0 y0 R", "cx y0 z0 R", "cy x0 z0 R").
            if mnemonic in ("cz", "cx", "cy"):
                axis = mnemonic[1]
                x0 = 0.0
                y0 = 0.0
                z0 = 0.0

                try:
                    if len(tokens) >= 5:
                        a0 = float(tokens[2])
                        b0 = float(tokens[3])
                        R = float(tokens[4])
                        if axis == "z":
                            x0, y0 = a0, b0
                        elif axis == "x":
                            y0, z0 = a0, b0
                        else:  # axis == "y"
                            x0, z0 = a0, b0
                    elif len(tokens) >= 3:
                        R = float(tokens[2])
                    else:
                        continue
                except ValueError:
                    continue

                if R <= 0.0:
                    continue

                dx = 1.0 if axis == "x" else 0.0
                dy = 1.0 if axis == "y" else 0.0
                dz = 1.0 if axis == "z" else 0.0
                cyl_list.append(
                    CylinderSurf(axis=axis, x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, R=R)
                )
                continue

            # Finite cylinder: RCC x0 y0 z0 dx dy dz R
            if mnemonic == "rcc" and len(tokens) >= 9:
                try:
                    x0 = float(tokens[2])
                    y0 = float(tokens[3])
                    z0 = float(tokens[4])
                    dx = float(tokens[5])
                    dy = float(tokens[6])
                    dz = float(tokens[7])
                    R = float(tokens[8])
                except ValueError:
                    continue
                if R <= 0.0:
                    continue
                axis = "z"
                if abs(dx) >= abs(dy) and abs(dx) >= abs(dz):
                    axis = "x"
                elif abs(dy) >= abs(dx) and abs(dy) >= abs(dz):
                    axis = "y"
                cyl_list.append(
                    CylinderSurf(axis=axis, x0=x0, y0=y0, z0=z0, dx=dx, dy=dy, dz=dz, R=R)
                )
                continue

    return px_list, py_list, pz_list, cyl_list, sph_list


def _build_infinite_box(gmsh, length=10.0):
    half = 0.5 * length
    tag = gmsh.model.occ.addBox(-half, -half, -half, length, length, length)
    gmsh.model.occ.synchronize()
    return [(3, tag)]


def _build_slab_box(gmsh, px_list, py_list, pz_list):
    def unique_sorted(values, tol=1.0e-12):
        out = []
        for v in sorted(values):
            if not out or abs(v - out[-1]) > tol:
                out.append(v)
        return out

    other_width = 1.0
    x_planes = unique_sorted([p.value for p in px_list])
    y_planes = unique_sorted([p.value for p in py_list])
    z_planes = unique_sorted([p.value for p in pz_list])

    axis = None
    cuts = []
    if len(x_planes) >= 2:
        axis = "x"
        cuts = x_planes
        y_min, y_max = -0.5 * other_width, 0.5 * other_width
        z_min, z_max = -0.5 * other_width, 0.5 * other_width
    elif len(y_planes) >= 2:
        axis = "y"
        cuts = y_planes
        x_min, x_max = -0.5 * other_width, 0.5 * other_width
        z_min, z_max = -0.5 * other_width, 0.5 * other_width
    elif len(z_planes) >= 2:
        axis = "z"
        cuts = z_planes
        x_min, x_max = -0.5 * other_width, 0.5 * other_width
        y_min, y_max = -0.5 * other_width, 0.5 * other_width
    else:
        # No axis has at least two planes; create a unit proxy box.
        tag = gmsh.model.occ.addBox(-0.5, -0.5, -0.5, 1.0, 1.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    regions = []
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        d = hi - lo
        if d <= 0.0:
            continue
        if axis == "x":
            tag = gmsh.model.occ.addBox(lo, y_min, z_min, d, y_max - y_min, z_max - z_min)
        elif axis == "y":
            tag = gmsh.model.occ.addBox(x_min, lo, z_min, x_max - x_min, d, z_max - z_min)
        else:
            tag = gmsh.model.occ.addBox(x_min, y_min, lo, x_max - x_min, y_max - y_min, d)
        regions.append((3, tag))

    if not regions:
        # Degenerate plane list after filtering; fallback box.
        tag = gmsh.model.occ.addBox(-0.5, -0.5, -0.5, 1.0, 1.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    if len(regions) > 1:
        # Force conformal interfaces between adjacent slab regions.
        # Without fragmentation, touching boxes can remain topologically
        # disconnected (coincident faces with distinct node sets).
        frag, _ = gmsh.model.occ.fragment(
            [regions[0]],
            regions[1:],
            removeObject=True,
            removeTool=True,
        )
        regions = [(dim, tag) for dim, tag in frag if dim == 3]
        gmsh.model.occ.synchronize()

        # Re-establish deterministic left-to-right (or y/z) region ordering.
        def axis_center(region):
            dim, tag = region
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
            if axis == "x":
                return 0.5 * (xmin + xmax)
            if axis == "y":
                return 0.5 * (ymin + ymax)
            return 0.5 * (zmin + zmax)

        regions = sorted(regions, key=axis_center)

    gmsh.model.occ.synchronize()
    return regions


def _build_sphere_regions(gmsh, spheres):
    if not spheres:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    radii = sorted({s.R for s in spheres if s.R > 0.0})
    if not radii:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, 1.0)
        gmsh.model.occ.synchronize()
        return [(3, tag)]

    sphere_tags = []
    for radius in radii:
        tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, radius)
        sphere_tags.append(tag)

    gmsh.model.occ.synchronize()

    regions = [(3, sphere_tags[0])]
    for i in range(1, len(sphere_tags)):
        outer = sphere_tags[i]
        inner = sphere_tags[i - 1]
        cut_vols, _ = gmsh.model.occ.cut(
            [(3, outer)],
            [(3, inner)],
            removeObject=True,
            removeTool=False,
        )
        for dim, tag in cut_vols:
            if dim == 3:
                regions.append((dim, tag))

    gmsh.model.occ.synchronize()
    return regions


def _build_cylinder_regions(gmsh, cyl_list, pz_list):
    if not cyl_list:
        return _build_slab_box(gmsh, [], [], [])

    c0 = cyl_list[0]
    axis = c0.axis

    if axis == "z" and pz_list:
        vals = sorted(p.value for p in pz_list)
        z_min, z_max = vals[0], vals[-1]
        base = (c0.x0, c0.y0, z_min)
        height = (0.0, 0.0, z_max - z_min)
    else:
        if abs(c0.dx) + abs(c0.dy) + abs(c0.dz) > 0.0:
            base = (c0.x0, c0.y0, c0.z0)
            height = (c0.dx, c0.dy, c0.dz)
        else:
            base = (c0.x0, c0.y0, -1.0)
            height = (0.0, 0.0, 2.0)

    radii = sorted({c.R for c in cyl_list if c.R > 0.0})
    cyl_tags = []
    for radius in radii:
        tag = gmsh.model.occ.addCylinder(
            base[0],
            base[1],
            base[2],
            height[0],
            height[1],
            height[2],
            radius,
        )
        cyl_tags.append(tag)

    gmsh.model.occ.synchronize()

    regions = [(3, cyl_tags[0])]
    for i in range(1, len(cyl_tags)):
        outer = cyl_tags[i]
        inner = cyl_tags[i - 1]
        cut_vols, _ = gmsh.model.occ.cut(
            [(3, outer)],
            [(3, inner)],
            removeObject=True,
            removeTool=False,
        )
        for dim, tag in cut_vols:
            if dim == 3:
                regions.append((dim, tag))

    gmsh.model.occ.synchronize()
    return regions


def _tag_axis_aligned_boundary_surfaces(gmsh):
    # Tag external axis-aligned faces so OpenSn boundary conditions can
    # be assigned deterministically.
    surfaces = gmsh.model.getEntities(2)
    if not surfaces:
        return

    bboxes = []
    for dim, tag in surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        bboxes.append((tag, xmin, ymin, zmin, xmax, ymax, zmax))

    gxmin = min(bb[1] for bb in bboxes)
    gymin = min(bb[2] for bb in bboxes)
    gzmin = min(bb[3] for bb in bboxes)
    gxmax = max(bb[4] for bb in bboxes)
    gymax = max(bb[5] for bb in bboxes)
    gzmax = max(bb[6] for bb in bboxes)

    xspan = gxmax - gxmin
    yspan = gymax - gymin
    zspan = gzmax - gzmin
    tol = 1.0e-7 * max(xspan, yspan, zspan, 1.0)

    face_tags = {
        "xmin": [],
        "xmax": [],
        "ymin": [],
        "ymax": [],
        "zmin": [],
        "zmax": [],
    }

    for tag, xmin, ymin, zmin, xmax, ymax, zmax in bboxes:
        x_flat = abs(xmax - xmin) <= tol
        y_flat = abs(ymax - ymin) <= tol
        z_flat = abs(zmax - zmin) <= tol

        if x_flat and abs(xmin - gxmin) <= tol:
            face_tags["xmin"].append(tag)
        elif x_flat and abs(xmax - gxmax) <= tol:
            face_tags["xmax"].append(tag)
        elif y_flat and abs(ymin - gymin) <= tol:
            face_tags["ymin"].append(tag)
        elif y_flat and abs(ymax - gymax) <= tol:
            face_tags["ymax"].append(tag)
        elif z_flat and abs(zmin - gzmin) <= tol:
            face_tags["zmin"].append(tag)
        elif z_flat and abs(zmax - gzmax) <= tol:
            face_tags["zmax"].append(tag)

    for name, tags in face_tags.items():
        if not tags:
            continue
        pg = gmsh.model.addPhysicalGroup(2, tags)
        gmsh.model.setPhysicalName(2, pg, name)


def mcnp_to_gmsh_mesh(mcnp_input, mesh_out="geom.msh"):
    try:
        import gmsh
    except ImportError as exc:
        raise SystemExit(
            "The 'gmsh' Python package is required for MCNP->Gmsh conversion.\n"
            "Install it with e.g. 'pip install gmsh'."
        ) from exc

    geom_kind = _detect_geometry_kind_from_name(mcnp_input)
    px_list, py_list, pz_list, cyl_list, sph_list = _parse_mcnp_surfaces(mcnp_input)

    gmsh.initialize()
    gmsh.model.add("mcnp_geometry")

    # Pick a geometry builder based on what the MCNP parser found (or what
    # the filename suffix hints the problem type should be).
    if sph_list:
        regions = _build_sphere_regions(gmsh, sph_list)
    elif cyl_list:
        regions = _build_cylinder_regions(gmsh, cyl_list, pz_list)
    elif px_list or py_list or pz_list or geom_kind in ("SL", "ISLC"):
        regions = _build_slab_box(gmsh, px_list, py_list, pz_list)
    elif geom_kind == "CY":
        regions = _build_cylinder_regions(gmsh, cyl_list, pz_list)
    elif geom_kind == "SP":
        regions = _build_sphere_regions(gmsh, sph_list)
    elif geom_kind == "IN":
        regions = _build_infinite_box(gmsh)
    else:
        regions = _build_slab_box(gmsh, px_list, py_list, pz_list)

    for i, region in enumerate(regions, start=1):
        dim, tag = region
        if dim != 3:
            continue
        gmsh.model.addPhysicalGroup(3, [tag], i)
        gmsh.model.setPhysicalName(3, i, f"region_{i}")

    _tag_axis_aligned_boundary_surfaces(gmsh)

    gmsh.model.mesh.generate(3)
    gmsh.write(mesh_out)
    gmsh.finalize()
    _normalize_text_file_to_lf(mesh_out)

    print(
        f"[MCNP->Gmsh] Geometry '{geom_kind}' -> {len(regions)} region(s), "
        f"mesh written to '{mesh_out}'."
    )


def make_mesh_for_problem(src_problem_dir, out_problem_dir, gmsh_exe="gmsh", msh_name="geom.msh"):
    # 1) Prefer explicit Gmsh geometry files.
    geo_files = [f for f in os.listdir(src_problem_dir) if f.lower().endswith(".geo")]
    step_files = [f for f in os.listdir(src_problem_dir) if f.lower().endswith(".step") or f.lower().endswith(".stp")]

    geom_src = None
    geom_type = None

    if geo_files:
        geom_src = os.path.join(src_problem_dir, geo_files[0])
        geom_type = "geo"
    elif step_files:
        geom_src = os.path.join(src_problem_dir, step_files[0])
        geom_type = "step"

    os.makedirs(out_problem_dir, exist_ok=True)
    msh_path = os.path.join(out_problem_dir, msh_name)

    if geom_src is not None:
        print(f" Geometry: using {geom_type.upper()} file {geom_src}")
        print(f" Geometry: Mesh -> {msh_path}")
        try:
            subprocess.check_call(
                [gmsh_exe, geom_src, "-3", "-o", msh_path],
                cwd=src_problem_dir,
            )
        except FileNotFoundError:
            raise SystemExit(
                f"Gmsh executable '{gmsh_exe}' not found. "
                "Install Gmsh or pass --gmsh with the full path."
            )
        _normalize_text_file_to_lf(msh_path)
        print(f" Geometry: created mesh {msh_path}\n")
        return

    # 2) Fallback: MCNP input -> Gmsh via Python API.
    # This route is more limited, but it lets the repo generate a usable mesh
    # even when only the MCNP input deck is available.
    mcnp_candidates = [
        f for f in os.listdir(src_problem_dir)
        if f.lower().endswith((".i", ".inp", ".mcnp"))
    ]
    if not mcnp_candidates:
        print(
            f" Geometry: no .geo/.step or MCNP input (.i/.inp/.mcnp) found in "
            f"{src_problem_dir}, skipping mesh.\n"
        )
        return

    mcnp_input = os.path.join(src_problem_dir, mcnp_candidates[0])
    print(f" Geometry: using MCNP input {mcnp_input}")
    print(f" Geometry: Mesh -> {msh_path}")

    mcnp_to_gmsh_mesh(mcnp_input, msh_path)
    print(f" Geometry: created mesh {msh_path}\n")
