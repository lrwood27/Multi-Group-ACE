import argparse
import os
import warnings

from .json_loader import load_from_json
from .converter import build_ace
from .ace_io import read_mg_ace, write_ace, write_summary, write_xsdir
from .hdf5_io import write_mgxs_hdf5_isotropic, write_problem_mgxs_library
from .gmsh_geometry import make_mesh_for_problem
from .continuous_energy import convert_ace_to_openmc_hdf5, detect_ace_type


def print_header_info(mg, nxs):
    print("\n=====> multi_group_ace: create simple multigroup ACE file\n")
    print(f"\t zaid = {mg.zaid}")
    print(f"\t za = {mg.za}")
    print(f"\t awr = {mg.awr}")
    print(f"\t tmp = {mg.tmp}")
    print(f"\t groups = {mg.NG}")
    print(f"\t xss size = {nxs[1]}")
    print(f"\t comment = {mg.comment_short}")
    print(f"\t date = {mg.date}\n")


def print_mcnp_snippet(mg, nxs):
    print("\n MCNP input snippet for this table:\n")
    print(f"mgopt f {mg.NG}")
    print(f"xs1 {mg.zaid} {mg.erg[0]:.5g} {mg.zaid} 0 1 1 {nxs[1]} 0 0 {mg.tmp:.6g}")
    print(f"m1 {mg.zaid} 1.0\n")


def convert_single_ace_to_hdf5(ace_file, h5_out=None, library_name=None):
    mg, _, _, _ = read_mg_ace(ace_file)
    if h5_out is None:
        base, _ = os.path.splitext(ace_file)
        h5_out = base + ".h5"
    write_mgxs_hdf5_isotropic(mg, h5_out, library_name=library_name)
    print(f"[HDF5] Wrote: {h5_out}")


def scan_verification_and_convert(root_dir, out_dir, library_name=None, make_mesh=False, gmsh_exe="gmsh"):
    root_dir = os.path.abspath(root_dir)
    out_dir = os.path.abspath(out_dir)
    if not os.path.isdir(root_dir):
        raise ValueError(f"Verification root is not a directory: {root_dir}")
    if os.path.normcase(root_dir) == os.path.normcase(out_dir):
        raise ValueError("The output directory must differ from the verification root.")
    if library_name is not None:
        raise ValueError("--library-name is supported only with --ace.")
    per_problem_mg = {}
    problem_sources = {}

    # Walk the tree, try to interpret files as ACE tables, and group any
    # successful reads by problem directory so each problem gets one mgxs.h5.
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = sorted(
            name for name in dirnames
            if os.path.normcase(os.path.join(dirpath, name)) != os.path.normcase(out_dir)
        )
        for fname in sorted(filenames):
            ace_path = os.path.join(dirpath, fname)
            try:
                mg, _, _, _ = read_mg_ace(ace_path)
            except (ValueError, OSError, UnicodeError) as exc:
                if fname.lower().endswith("m") and "." in fname:
                    warnings.warn(f"Skipping invalid ACE file {ace_path}: {exc}")
                continue

            rel = os.path.relpath(ace_path, root_dir)
            rel_dir = os.path.dirname(rel)
            out_problem_dir = os.path.join(out_dir, rel_dir)
            os.makedirs(out_problem_dir, exist_ok=True)

            key = os.path.normpath(out_problem_dir)
            per_problem_mg.setdefault(key, []).append(mg)
            problem_sources[key] = dirpath

    if not per_problem_mg:
        raise ValueError(f"No supported multigroup ACE tables found under {root_dir}.")

    for out_problem_dir, mg_list in per_problem_mg.items():
        write_problem_mgxs_library(mg_list, os.path.join(out_problem_dir, "mgxs.h5"))
        if make_mesh:
            make_mesh_for_problem(problem_sources[out_problem_dir], out_problem_dir, gmsh_exe=gmsh_exe)

    print(f"Scan complete. Wrote {len(per_problem_mg)} mgxs.h5 file(s).")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Multigroup ACE utilities and continuous-energy neutron ACE to OpenMC HDF5 conversion")

    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--config", "-config", "-j", help="JSON config for JSON -> multigroup ACE")
    modes.add_argument("--ace", nargs="+", metavar="FILE", help="ACE input(s): one multigroup table, or continuous-energy neutron libraries")
    parser.add_argument("--h5", help="Output HDF5 path")
    parser.add_argument("--library-name", default=None, help="Optional multigroup HDF5 library group name")

    modes.add_argument("--verification-root", help="Root directory for recursive multigroup ACE scan")
    parser.add_argument("--outdir", default=None, help="Output root for batch MGXS or continuous-energy libraries (default: h5_out)")
    parser.add_argument("--make-mesh", action="store_true", help="Also generate geometry mesh")
    parser.add_argument("--gmsh", default="gmsh", help="Gmsh executable path")
    parser.add_argument("--metastable-scheme", choices=("mcnp", "nndc"), default=None, help="Continuous-energy ZAID convention (default: mcnp)")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing continuous-energy HDF5 files and XML catalog")

    args = parser.parse_args(argv)
    if (args.h5 or args.library_name is not None) and not args.ace:
        parser.error("--h5 and --library-name require --ace.")
    if args.make_mesh and not args.verification_root:
        parser.error("--make-mesh requires --verification-root.")
    if (args.overwrite or args.metastable_scheme is not None) and not args.ace:
        parser.error("--overwrite and --metastable-scheme require continuous-energy --ace inputs.")
    if args.outdir is not None and args.config:
        parser.error("--outdir requires --ace or --verification-root.")

    # The CLI supports three modes:
    #   verification-root -> many ACE tables to many mgxs.h5 files (+ mesh)
    #   ace               -> MGXS or continuous-energy neutron HDF5 + XML
    #   config            -> JSON description to ACE output
    if args.verification_root:
        scan_verification_and_convert(
            root_dir=args.verification_root,
            out_dir=args.outdir or "h5_out",
            library_name=args.library_name,
            make_mesh=args.make_mesh,
            gmsh_exe=args.gmsh,
        )
        return

    if args.ace:
        try:
            kinds = {detect_ace_type(path) for path in args.ace}
            if len(kinds) != 1:
                raise ValueError("Cannot mix multigroup, continuous-energy, or unsupported ACE table types in one conversion.")
            kind = kinds.pop()
            if kind == "multigroup":
                if len(args.ace) != 1:
                    raise ValueError("Multigroup --ace accepts one input; use --verification-root for batch conversion.")
                if args.overwrite or args.metastable_scheme is not None:
                    raise ValueError("--overwrite and --metastable-scheme are continuous-energy-only options.")
                if args.outdir is not None:
                    raise ValueError("Use --h5 for single-file multigroup conversion; --outdir is for libraries.")
                convert_single_ace_to_hdf5(args.ace[0], h5_out=args.h5, library_name=args.library_name)
            elif kind == "continuous-energy":
                if args.library_name is not None:
                    raise ValueError("--library-name is multigroup-only; OpenMC determines continuous-energy nuclide names.")
                if args.h5 and args.outdir is not None:
                    raise ValueError("--h5 cannot be combined with an explicit --outdir.")
                paths, xml = convert_ace_to_openmc_hdf5(
                    args.ace, args.outdir or "h5_out", h5_path=args.h5,
                    metastable_scheme=args.metastable_scheme or "mcnp",
                    overwrite=args.overwrite,
                )
                for path in [*paths, xml]:
                    print(f"[OpenMC] Wrote: {path}")
            else:
                raise ValueError("Unsupported ACE table type: use multigroup 'm' or continuous-energy neutron 'c' tables. Thermal scattering is not supported.")
        except (ValueError, OSError, ImportError, NotImplementedError) as exc:
            parser.error(str(exc))
        return

    if not args.config:
        parser.error("Provide --config, or --ace, or --verification-root.")

    # JSON-to-ACE path: normalize JSON input, build ACE arrays, then emit
    # the table plus the helper files used for inspection/MCNP setup.
    mg, ace_file = load_from_json(args.config)
    nxs, jxs, xss = build_ace(mg)

    print_header_info(mg, nxs)
    write_summary(ace_file, mg, nxs, jxs, xss)
    write_xsdir("xsdir_mcnp6.3", mg, ace_file, nxs)
    write_ace(ace_file, mg, nxs, jxs, xss)
    print_mcnp_snippet(mg, nxs)


if __name__ == "__main__":
    main()
