import argparse
import os

from .json_loader import load_from_json
from .converter import build_ace
from .ace_io import read_mg_ace, write_ace, write_summary, write_xsdir
from .hdf5_io import write_mgxs_hdf5_isotropic, write_problem_mgxs_library
from .gmsh_geometry import make_mesh_for_problem


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
    per_problem_mg = {}
    mesh_done_dirs = set()

    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            ace_path = os.path.join(dirpath, fname)
            try:
                mg, _, _, _ = read_mg_ace(ace_path)
            except Exception:
                continue

            rel = os.path.relpath(ace_path, root_dir)
            rel_dir = os.path.dirname(rel)
            out_problem_dir = os.path.join(out_dir, rel_dir)
            os.makedirs(out_problem_dir, exist_ok=True)

            key = os.path.normpath(out_problem_dir)
            per_problem_mg.setdefault(key, []).append(mg)

            if make_mesh:
                src_key = os.path.normpath(dirpath)
                if src_key not in mesh_done_dirs:
                    make_mesh_for_problem(dirpath, out_problem_dir, gmsh_exe=gmsh_exe)
                    mesh_done_dirs.add(src_key)

    for out_problem_dir, mg_list in per_problem_mg.items():
        write_problem_mgxs_library(mg_list, os.path.join(out_problem_dir, "mgxs.h5"))

    print(f"Scan complete. Wrote {len(per_problem_mg)} mgxs.h5 file(s).")


def main():
    parser = argparse.ArgumentParser(description="Multi-Group ACE processor")

    parser.add_argument("--config", help="JSON config for JSON -> ACE")
    parser.add_argument("--ace", help="ACE input for ACE -> HDF5")
    parser.add_argument("--h5", help="Output HDF5 path")
    parser.add_argument("--library-name", default=None, help="Optional HDF5 library group name")

    parser.add_argument("--verification-root", help="Root directory for recursive ACE scan")
    parser.add_argument("--outdir", default="h5_out", help="Output root for scanned HDF5 files")
    parser.add_argument("--make-mesh", action="store_true", help="Also generate geometry mesh")
    parser.add_argument("--gmsh", default="gmsh", help="Gmsh executable path")

    args = parser.parse_args()

    if args.verification_root:
        scan_verification_and_convert(
            root_dir=args.verification_root,
            out_dir=args.outdir,
            library_name=args.library_name,
            make_mesh=args.make_mesh,
            gmsh_exe=args.gmsh,
        )
        return

    if args.ace:
        convert_single_ace_to_hdf5(args.ace, h5_out=args.h5, library_name=args.library_name)
        return

    if not args.config:
        parser.error("Provide --config, or --ace, or --verification-root.")

    mg, ace_file = load_from_json(args.config)
    nxs, jxs, xss = build_ace(mg)

    print_header_info(mg, nxs)
    write_summary(ace_file, mg, nxs, jxs, xss)
    write_xsdir("xsdir_mcnp6.3", mg, ace_file, nxs)
    write_ace(ace_file, mg, nxs, jxs, xss)
    print_mcnp_snippet(mg, nxs)


if __name__ == "__main__":
    main()
