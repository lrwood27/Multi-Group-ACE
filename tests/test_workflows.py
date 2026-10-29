"""Regression checks for the maintained multigroup conversion workflows.

Run from the repository root: python -m unittest discover -s tests -v
"""

import contextlib
import copy
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from multi_group_ace import build_ace, load_from_json, read_mg_ace, write_ace
from multi_group_ace.ace_io import _read_control_array
from multi_group_ace.cli import scan_verification_and_convert
from multi_group_ace.converter import get_equiprob_bins
from multi_group_ace.core import KB_EV_PER_K, MEV_TO_EV
from multi_group_ace.gmsh_geometry import make_mesh_for_problem
from multi_group_ace.hdf5_io import write_mgxs_hdf5_isotropic, write_problem_mgxs_library


ROOT = Path(__file__).resolve().parents[1]
ACE_FILES = sorted(p for p in (ROOT / "examples").glob("mg*/*") if p.suffix.endswith("m"))
HAS_H5PY = importlib.util.find_spec("h5py") is not None
HAS_GMSH = importlib.util.find_spec("gmsh") is not None


class WorkflowCase(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def material(self):
        return load_from_json(ROOT / "json examples" / "ex4.json")[0]

    def run_cli(self, *args, legacy=False):
        entry = [str(ROOT / "simple_ace_mg_full.py")] if legacy else ["-m", "multi_group_ace.cli"]
        env = dict(os.environ, PYTHONPATH=str(ROOT))
        return subprocess.run(
            [sys.executable, *entry, *map(str, args)], cwd=self.work,
            env=env, text=True, capture_output=True, timeout=30,
        )


class ACEChecks(WorkflowCase):
    def test_bundled_tables_round_trip_and_conserve_neutrons(self):
        self.assertEqual(len(ACE_FILES), 15)
        for source in ACE_FILES:
            with self.subTest(source=source):
                mg = read_mg_ace(source)[0]
                for incoming in range(mg.NG):
                    scatter = sum(mg.sigs[outgoing * mg.NG + incoming] for outgoing in range(mg.NG))
                    self.assertAlmostEqual(mg.sigt[incoming], mg.sigc[incoming] + mg.sigf[incoming] + scatter)
                arrays = build_ace(mg)
                target = self.work / "roundtrip.ace"
                write_ace(target, mg, *arrays)
                restored = read_mg_ace(target)[0]
                for field in ("erg", "sigt", "sigc", "sigf", "nu", "chi", "sigs", "sig1"):
                    for expected, actual in zip(getattr(mg, field), getattr(restored, field)):
                        self.assertAlmostEqual(expected, actual, delta=1e-9)

    @unittest.skipUnless(importlib.util.find_spec("numpy"), "numpy is required for benchmark checks")
    def test_infinite_medium_benchmarks(self):
        # Independent zero-leakage balance: (diag(total) - scatter) phi = F phi/k.
        import numpy as np
        for name in ("mg01", "mg05", "mg31", "mg74", "mg75"):
            with self.subTest(problem=name):
                source = next(p for p in ACE_FILES if p.parent.name == name)
                mg = read_mg_ace(source)[0]
                scatter = np.array(mg.sigs).reshape(mg.NG, mg.NG)
                loss = np.diag(mg.sigt) - scatter
                production = np.outer(mg.chi, np.array(mg.nu) * mg.sigf)
                k = max(np.linalg.eigvals(np.linalg.solve(loss, production)).real)
                description = json.loads((source.parent / "description.json").read_text())
                self.assertAlmostEqual(k, description["analytic_data"]["k-eff"]["val"], delta=1e-6)

    def test_empty_input_fails(self):
        target = self.work / "empty"
        target.touch()
        with self.assertRaisesRegex(ValueError, "multigroup ACE"):
            read_mg_ace(target)

    def test_truncated_control_arrays_fail_at_eof(self):
        for count, name in ((16, "NXS"), (32, "JXS")):
            with self.subTest(array=name), self.assertRaisesRegex(ValueError, name):
                _read_control_array(io.StringIO("0 1\n\n"), count, name)

    def test_continuous_energy_input_is_rejected_explicitly(self):
        target = self.work / "U235.ace"
        target.write_text("92235.80c 233.0248 2.5301e-8\ncomment\n")
        with self.assertRaisesRegex(ValueError, "Continuous-energy ACE tables are not supported"):
            read_mg_ace(target)

    def test_truncated_payload_and_invalid_pointers_fail(self):
        mg = self.material()
        arrays = build_ace(mg)
        target = self.work / "bad.ace"
        write_ace(target, mg, *arrays)
        target.write_text("\n".join(target.read_text().splitlines()[:-1]))
        with self.assertRaisesRegex(ValueError, "XSS"):
            read_mg_ace(target)
        nxs, jxs, xss = build_ace(mg)
        jxs[2] = nxs[1] + 1
        write_ace(target, mg, nxs, jxs, xss)
        with self.assertRaisesRegex(ValueError, "outside"):
            read_mg_ace(target)

    def test_anisotropic_round_trip(self):
        mg = self.material()
        mg.sig1 = [0.08, 0.0, -0.04, 0.04]
        target = self.work / "p1.ace"
        write_ace(target, mg, *build_ace(mg))
        restored = read_mg_ace(target)[0]
        for expected, actual in zip(mg.sig1, restored.sig1):
            self.assertAlmostEqual(expected, actual, delta=1e-5)


class AngularChecks(unittest.TestCase):
    def test_bins_invert_the_analytic_linear_cdf(self):
        for nbins in (10, 100, 1000, 2000):
            for moment in (-0.2, 0.0, 0.2):
                with self.subTest(nbins=nbins, moment=moment):
                    mu = get_equiprob_bins(nbins, 1, [1.0, moment])
                    self.assertEqual(mu[0], -1.0)
                    self.assertEqual(mu[-1], 1.0)
                    self.assertTrue(all(a < b for a, b in zip(mu, mu[1:])))
                    for k, x in enumerate(mu):
                        cdf = 0.5 * (1 + x) + 0.75 * moment * (x*x - 1)
                        self.assertAlmostEqual(cdf, k/nbins, delta=1e-6)

    def test_nonpositive_bin_count_fails(self):
        for nbins in (0, -1):
            with self.subTest(nbins=nbins), self.assertRaisesRegex(ValueError, "positive integer"):
                get_equiprob_bins(nbins, 1, [1.0, 0.1])

    def test_negative_angular_probability_fails(self):
        with self.assertRaisesRegex(ValueError, "negative angular"):
            get_equiprob_bins(100, 1, [1.0, 0.4])


@unittest.skipUnless(HAS_H5PY, "h5py is required for HDF5 integration checks")
class HDF5Checks(WorkflowCase):
    def test_json_cli_to_ace_to_hdf5(self):
        import h5py
        result = self.run_cli("--config", ROOT / "json examples" / "ex4.json")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.work / "xsdir_mcnp6.3").is_file())
        self.assertTrue((self.work / "mg02xsecs.summary.txt").is_file())
        result = self.run_cli("--ace", self.work / "mg02xsecs", "--h5", self.work / "nested" / "mgxs.h5", "--library-name", "material")
        self.assertEqual(result.returncode, 0, result.stderr)
        with h5py.File(self.work / "nested" / "mgxs.h5") as f:
            self.assertEqual(f.attrs["filetype"], b"mgxs")
            self.assertEqual(f.attrs["energy_groups"], 2)
            self.assertEqual(f["material/294K/scatter_data/scatter_matrix"][:].tolist(), [0.4, 0.4, 0.0, 0.4])
            for actual, expected in zip(f.attrs["group structure"], [10.0, 625000.0, 20000000.0]):
                self.assertAlmostEqual(actual, expected, delta=1e-8)
            self.assertEqual(f["material/294K/absorption"][:].tolist(), [0.2, 0.1])
            self.assertNotIn("fission", f["material/294K"])

    def test_p1_metadata_and_scattering_layout(self):
        import h5py
        mg = self.material()
        mg.sig1 = [0.08, 0.0, -0.04, 0.04]
        target = self.work / "p1.h5"
        write_mgxs_hdf5_isotropic(mg, target)
        with h5py.File(target) as f:
            group = f[mg.zaid]
            self.assertEqual(group.attrs["representation"], b"isotropic")
            self.assertEqual(group.attrs["order"], 1)
            self.assertEqual(group.attrs["scatter_format"], b"legendre")
            self.assertEqual(group.attrs["scatter_shape"], b"[G][G'][Order]")
            self.assertEqual(group["294K/scatter_data/scatter_matrix"][:].tolist(), [0.4, 0.08, 0.4, -0.04, 0., 0., 0.4, 0.04])

    def test_materials_keep_individual_temperatures(self):
        import h5py
        mg = self.material()
        other = copy.deepcopy(mg)
        other.zaid = "99954.01m"
        other.tmp = 600 * KB_EV_PER_K / MEV_TO_EV
        target = self.work / "combined.h5"
        write_problem_mgxs_library([mg, other], target)
        with h5py.File(target) as f:
            self.assertIn("294K", f[mg.zaid])
            self.assertIn("600K", f[other.zaid])
            self.assertAlmostEqual(f[other.zaid]["kTs/600K"][()], 600 * KB_EV_PER_K)

    def test_zero_temperature_has_zero_kelvin_label(self):
        import h5py
        mg = self.material()
        mg.tmp = 0.
        target = self.work / "zero.h5"
        write_mgxs_hdf5_isotropic(mg, target)
        with h5py.File(target) as f:
            self.assertEqual(f[mg.zaid]["kTs/0K"][()], 0.)

    def test_incompatible_boundaries_are_rejected_before_writing(self):
        mg = self.material()
        other = copy.deepcopy(mg)
        other.zaid = "99954.01m"
        other.erg[1] = 1.0
        target = self.work / "combined.h5"
        target.write_bytes(b"existing output")
        with self.assertRaisesRegex(ValueError, "Energy grids differ"):
            write_problem_mgxs_library([mg, other], target)
        self.assertEqual(target.read_bytes(), b"existing output")

    def test_duplicate_names_are_rejected_before_writing(self):
        mg = self.material()
        target = self.work / "duplicate.h5"
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            write_problem_mgxs_library([mg, mg], target)
        self.assertFalse(target.exists())

    def test_batch_conversion_of_all_bundled_tables(self):
        import h5py
        target = self.work / "verification"
        with contextlib.redirect_stdout(io.StringIO()):
            scan_verification_and_convert(ROOT / "examples", target)
        outputs = sorted(target.glob("*/mgxs.h5"))
        self.assertEqual(len(outputs), 13)
        count = 0
        for output in outputs:
            with h5py.File(output) as f:
                sources = [p for p in ACE_FILES if p.parent.name == output.parent.name]
                self.assertEqual(sorted(f), sorted(read_mg_ace(p)[0].zaid for p in sources))
                for source in sources:
                    mg = read_mg_ace(source)[0]
                    group = f[mg.zaid]
                    temperature = next(iter(group["kTs"]))
                    data = group[temperature]
                    self.assertEqual(data["total"][:].tolist(), mg.sigt)
                    self.assertEqual(data["absorption"][:].tolist(), [c+sf for c,sf in zip(mg.sigc,mg.sigf)])
                    if any(mg.sigf):
                        self.assertEqual(data["nu-fission"][:].tolist(), [n*sf for n,sf in zip(mg.nu,mg.sigf)])
                    count += 1
        self.assertEqual(count, 15)

    def test_batch_skips_its_own_output_tree(self):
        root = self.work / "sources"
        problem = root / "problem"
        problem.mkdir(parents=True)
        (problem / "table.01m").write_bytes(ACE_FILES[0].read_bytes())
        target = root / "output"
        (target / "stale").mkdir(parents=True)
        (target / "stale" / "table.01m").write_bytes(ACE_FILES[0].read_bytes())
        with contextlib.redirect_stdout(io.StringIO()):
            scan_verification_and_convert(root, target)
        self.assertTrue((target / "problem" / "mgxs.h5").exists())
        self.assertFalse((target / "output").exists())

    def test_legacy_cli_uses_maintained_conversion(self):
        result = self.run_cli("--config", ROOT / "json examples" / "ex4.json", legacy=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.work / "mg02xsecs").is_file())

    @unittest.skipUnless(HAS_GMSH, "gmsh is required for mesh integration checks")
    def test_batch_cli_generates_all_bundled_meshes(self):
        import gmsh
        target = self.work / "verification"
        result = self.run_cli("--verification-root", ROOT / "examples", "--outdir", target, "--make-mesh")
        self.assertEqual(result.returncode, 0, result.stderr)
        meshes = sorted(target.glob("*/geom.msh"))
        self.assertEqual(len(meshes), 13)
        self.assertEqual(len(list(target.glob("*/mgxs.h5"))), 13)
        for mesh in meshes:
            with self.subTest(mesh=mesh):
                gmsh.initialize()
                try:
                    gmsh.option.setNumber("General.Terminal", 0)
                    gmsh.open(str(mesh))
                    self.assertGreater(len(gmsh.model.mesh.getNodes()[0]), 0)
                    self.assertGreater(sum(len(tags) for tags in gmsh.model.mesh.getElements(3)[1]), 0)
                    self.assertGreater(len(gmsh.model.getPhysicalGroups(3)), 0)
                finally:
                    gmsh.finalize()


class InterfaceChecks(WorkflowCase):
    def test_conflicting_cli_modes_are_rejected(self):
        result = self.run_cli("--config", "input.json", "--ace", "input.ace")
        self.assertEqual(result.returncode, 2)
        self.assertIn("not allowed", result.stderr)

    def test_single_file_options_are_rejected_in_batch_mode(self):
        result = self.run_cli("--verification-root", ROOT / "examples", "--library-name", "ignored")
        self.assertEqual(result.returncode, 2)

    def test_invalid_or_empty_scan_root_is_rejected(self):
        for source in (self.work / "missing", self.work):
            with self.subTest(source=source), self.assertRaises(ValueError):
                scan_verification_and_convert(source, self.work / "output")

    def test_explicit_geometry_uses_absolute_paths(self):
        source = self.work / "source"
        source.mkdir()
        (source / "cube.geo").write_text('SetFactory("OpenCASCADE"); Box(1) = {0,0,0,1,1,1};')
        target = self.work / "output"

        def fake_gmsh(command, cwd):
            self.assertEqual(command[0], "resolved-gmsh.bat")
            self.assertTrue(Path(command[1]).is_absolute())
            self.assertTrue(Path(command[-1]).is_absolute())
            self.assertEqual(Path(cwd), source)
            Path(command[-1]).write_bytes(b"$MeshFormat\r\n2.2 0 8\r\n")

        with patch("multi_group_ace.gmsh_geometry.shutil.which", return_value="resolved-gmsh.bat"), patch(
            "multi_group_ace.gmsh_geometry.subprocess.check_call", side_effect=fake_gmsh
        ):
            with contextlib.redirect_stdout(io.StringIO()):
                make_mesh_for_problem(os.path.relpath(source), os.path.relpath(target))
        self.assertNotIn(b"\r", (target / "geom.msh").read_bytes())


if __name__ == "__main__":
    unittest.main()
