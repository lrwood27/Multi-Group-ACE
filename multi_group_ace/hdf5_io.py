from .core import KB_EV_PER_K, MEV_TO_EV
from math import isclose, isfinite
from pathlib import Path


def _set_fixed_ascii_attr(attrs, key, value):
    # OpenSn expects fixed-length C-string attributes.
    encoded = value.encode("ascii")
    attrs.create(key, encoded, dtype=f"S{len(encoded)}")


def _temp_group_name(tmp_mev):
    kT_eV = tmp_mev * MEV_TO_EV
    if not isfinite(kT_eV) or kT_eV < 0.0:
        raise ValueError("Temperature must be finite and nonnegative.")
    T_K = kT_eV / KB_EV_PER_K
    return f"{int(round(T_K))}K", kT_eV


def _scatter_order(mg):
    sig1 = getattr(mg, "sig1", [])
    if any(abs(x) > 0.0 for x in sig1):
        return 1
    return 0


def _scatter_matrix_openmc_order(mg, order):
    # Internal matrix convention:
    #   mg.sigs[row, col] = Sigma_s(from=col -> to=row), row-major.
    # OpenMC/OpenSn flattening convention for "scatter_matrix":
    #   incoming group (g_in) outer, outgoing group (g_out) inner, then order.
    # For P1, each (g_in, g_out) contributes [P0, P1].
    ordered = []
    for g_in in range(mg.NG):
        for g_out in range(mg.NG):
            idx = g_out * mg.NG + g_in
            ordered.append(mg.sigs[idx])
            if order >= 1:
                ordered.append(mg.sig1[idx])
    return ordered


def _write_nuclide_group(root, name, mg, temp_group_name, kT_eV):
    fissionable = any(sf > 0.0 for sf in mg.sigf)
    order = _scatter_order(mg)
    # This describes dependence on the incoming direction, independently of
    # the Legendre expansion of the outgoing scattering distribution.
    representation = "isotropic"

    # Each material/nuclide gets its own group under the file root. The
    # temperature-specific data then live inside a "<TTT>K" subgroup.
    g_lib = root.create_group(name)
    g_lib.attrs["atomic_weight_ratio"] = mg.awr
    g_lib.attrs["fissionable"] = bool(fissionable)
    _set_fixed_ascii_attr(g_lib.attrs, "representation", representation)
    _set_fixed_ascii_attr(g_lib.attrs, "scatter_format", "legendre")
    _set_fixed_ascii_attr(g_lib.attrs, "scatter_shape", "[G][G'][Order]")
    g_lib.attrs["order"] = order

    g_kTs = g_lib.create_group("kTs")
    g_kTs.create_dataset(temp_group_name, data=kT_eV)

    g_T = g_lib.create_group(temp_group_name)
    g_T.create_dataset("total", data=mg.sigt)
    g_T.create_dataset("absorption", data=[c + f for c, f in zip(mg.sigc, mg.sigf)])

    if fissionable:
        g_T.create_dataset("fission", data=mg.sigf)
        g_T.create_dataset("chi", data=mg.chi)
        g_T.create_dataset("nu-fission", data=[nu * f for nu, f in zip(mg.nu, mg.sigf)])

    # OpenMC/OpenSn-style scatter data are stored in a dedicated subgroup.
    g_scat = g_T.create_group("scatter_data")
    g_scat.create_dataset("g_min", data=[1] * mg.NG)
    g_scat.create_dataset("g_max", data=[mg.NG] * mg.NG)
    g_scat.create_dataset("scatter_matrix", data=_scatter_matrix_openmc_order(mg, order))
    g_scat.create_dataset("multiplicity_matrix", data=[1.0] * (mg.NG * mg.NG))


def write_mgxs_hdf5_isotropic(mg, h5_path, library_name=None):
    try:
        import h5py
    except ImportError as exc:
        raise SystemExit("h5py is required. Install with: pip install h5py") from exc

    temp_group_name, kT_eV = _temp_group_name(mg.tmp)
    # The HDF5 schema expects the group structure in ascending eV, while the
    # package stores energies internally as descending MeV boundaries.
    e_asc_eV = [e * MEV_TO_EV for e in reversed(mg.erg)]
    Path(h5_path).parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "w") as f:
        # Root attributes describe the library as a whole.
        _set_fixed_ascii_attr(f.attrs, "filetype", "mgxs")
        f.attrs["version"] = [1, 0]
        f.attrs["energy_groups"] = mg.NG
        f.attrs["delayed_groups"] = 0
        f.attrs["group structure"] = e_asc_eV

        _write_nuclide_group(f, library_name or mg.zaid, mg, temp_group_name, kT_eV)


def write_problem_mgxs_library(mg_list, h5_path):
    if not mg_list:
        return

    try:
        import h5py
    except ImportError as exc:
        raise SystemExit("h5py is required. Install with: pip install h5py") from exc

    ref = mg_list[0]
    for mg in mg_list[1:]:
        if mg.NG != ref.NG:
            raise ValueError("All ACE files in one problem must have same NG.")
        if len(mg.erg) != len(ref.erg) or not all(
            isclose(a, b, rel_tol=1.0e-10, abs_tol=1.0e-12)
            for a, b in zip(mg.erg, ref.erg)
        ):
            raise ValueError("Energy grids differ between ACE files.")

    names = [mg.zaid for mg in mg_list]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate ACE table names in one problem library.")
    temperatures = [_temp_group_name(mg.tmp) for mg in mg_list]

    # All materials in one problem library must share the same global energy
    # grid so a single mesh can reference them consistently.
    e_asc_eV = [e * MEV_TO_EV for e in reversed(ref.erg)]
    Path(h5_path).parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "w") as f:
        _set_fixed_ascii_attr(f.attrs, "filetype", "mgxs")
        f.attrs["version"] = [1, 0]
        f.attrs["energy_groups"] = ref.NG
        f.attrs["delayed_groups"] = 0
        f.attrs["group structure"] = e_asc_eV

        for mg, (temp_group_name, kT_eV) in zip(mg_list, temperatures):
            _write_nuclide_group(f, mg.zaid, mg, temp_group_name, kT_eV)
