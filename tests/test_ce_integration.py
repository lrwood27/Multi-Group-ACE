"""Real OpenMC tests using checksum-pinned authentic ACE data.

Prepare with ``python tests/fetch_ce_data.py`` and install requirements-ce.txt.
Set MULTI_GROUP_ACE_REQUIRE_CE_TESTS=1 to turn missing prerequisites into errors.
Synthetic temperature tests change only ACE metadata; they are plumbing checks,
not independently processed physical cross sections at the synthetic temperature.
"""

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from multi_group_ace.continuous_energy import convert_ace_to_openmc_hdf5


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / ".test-data" / "ce"
MANIFEST = json.loads(Path(__file__).with_name("ce_data_manifest.json").read_text())


class ContinuousEnergyIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        missing = []
        if importlib.util.find_spec("openmc") is None:
            missing.append("OpenMC Python API (install requirements-ce.txt)")
        for entry in MANIFEST["files"]:
            path = DATA / entry["name"]
            if not path.is_file():
                missing.append(f"fixture {path} (run python tests/fetch_ce_data.py)")
            elif path.stat().st_size != entry["size"] or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
                raise RuntimeError(f"Authentic fixture failed checksum verification: {path}")
        if missing:
            message = "; ".join(missing)
            if os.environ.get("MULTI_GROUP_ACE_REQUIRE_CE_TESTS") == "1":
                raise RuntimeError(message)
            raise unittest.SkipTest(message)
        import openmc.data
        cls.data_api = openmc.data

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.he3 = DATA / "2003.710nc"
        self.u235 = DATA / "92235.710nc"

    def compare_hdf5(self, actual, expected):
        """Compare every persisted dataset and attribute, including distributions."""
        import h5py
        import numpy as np
        with h5py.File(actual) as left, h5py.File(expected) as right:
            def compare(a, b, name):
                self.assertEqual(set(a.attrs), set(b.attrs), name)
                for key in a.attrs:
                    np.testing.assert_array_equal(a.attrs[key], b.attrs[key], err_msg=f"{name}@{key}")
                if isinstance(a, h5py.Dataset):
                    np.testing.assert_array_equal(a[()], b[()], err_msg=name)
                else:
                    self.assertEqual(set(a), set(b), name)
                    for key in a:
                        compare(a[key], b[key], f"{name}/{key}")
            compare(left, right, "")

    def assert_catalog(self, xml, expected_names):
        catalog = self.data_api.DataLibrary.from_xml(xml)
        self.assertEqual(len(catalog), len(expected_names))
        for name in expected_names:
            entry = catalog.get_by_material(name)
            self.assertIsNotNone(entry, name)
            self.assertEqual(entry["type"], "neutron")
            self.assertTrue(Path(entry["path"]).is_file(), entry)

    def synthetic_he3_temperature(self, kelvin=600):
        """Change only metadata to exercise merging, not Doppler broadening."""
        lines = self.he3.read_text().splitlines(keepends=True)
        kT_mev = kelvin * self.data_api.K_BOLTZMANN / 1e6
        for index in (1, 3):
            lines[index] = lines[index].replace("2.5301E-08", f"{kT_mev:.4E}")
        target = self.work / f"synthetic-He3-{kelvin}K.ace"
        target.write_text("".join(lines))
        return target

    def test_authentic_he3_and_u235_cli_round_trip(self):
        import h5py
        import numpy as np
        target = self.work / "converted"
        result = subprocess.run(
            [sys.executable, "-m", "multi_group_ace.cli", "--ace", str(self.he3), str(self.u235), "--outdir", str(target)],
            cwd=ROOT, capture_output=True, text=True, timeout=120,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for source in (self.he3, self.u235):
            direct = self.data_api.IncidentNeutron.from_ace(str(source), metastable_scheme="mcnp")
            actual = target / f"{direct.name}.h5"
            restored = self.data_api.IncidentNeutron.from_hdf5(actual)
            self.assertEqual(restored.name, direct.name)
            self.assertEqual(restored.temperatures, direct.temperatures)
            self.assertEqual(restored.atomic_number, direct.atomic_number)
            self.assertEqual(restored.mass_number, direct.mass_number)
            np.testing.assert_array_equal(restored.kTs, direct.kTs)
            for temperature in direct.temperatures:
                np.testing.assert_array_equal(restored.energy[temperature], direct.energy[temperature])
            reference = self.work / f"reference-{direct.name}.h5"
            direct.export_to_hdf5(reference, mode="w")
            self.compare_hdf5(actual, reference)
            with h5py.File(actual) as f:
                self.assertEqual(f.attrs["filetype"], b"data_neutron")
                self.assertIn("reactions", f[direct.name])
            if direct.name == "U235":
                self.assertIn(18, restored.reactions)
                self.assertGreater(len(restored.reactions[18].products), 0)
        self.assert_catalog(target / "cross_sections.xml", {"He3", "U235"})

    def test_all_tables_in_a_concatenated_library_are_converted(self):
        combined = self.work / "library.without_ace_suffix"
        combined.write_bytes(self.he3.read_bytes().rstrip() + b"\n" + self.u235.read_bytes())
        paths, xml = convert_ace_to_openmc_hdf5(combined, self.work / "out")
        self.assertEqual({p.name for p in paths}, {"He3.h5", "U235.h5"})
        self.assert_catalog(xml, {"He3", "U235"})

    def test_binary_ace_uses_the_real_openmc_parser(self):
        binary = self.work / "binary.dat"
        self.data_api.ace.ascii_to_binary(self.he3, binary)
        paths, xml = convert_ace_to_openmc_hdf5(binary, self.work / "binary-out")
        reference = self.work / "reference.h5"
        self.data_api.IncidentNeutron.from_ace(self.he3, metastable_scheme="mcnp").export_to_hdf5(reference, mode="w")
        self.compare_hdf5(paths[0], reference)
        self.assert_catalog(xml, {"He3"})

    def test_synthetic_temperature_metadata_merges_and_retains_cross_sections(self):
        import numpy as np
        second = self.synthetic_he3_temperature()
        paths, xml = convert_ace_to_openmc_hdf5([second, self.he3], self.work / "out")
        self.assertEqual(len(paths), 1)
        restored = self.data_api.IncidentNeutron.from_hdf5(paths[0])
        self.assertEqual(restored.temperatures, ["294K", "600K"])
        np.testing.assert_array_equal(restored.energy["294K"], restored.energy["600K"])
        for reaction in restored.reactions.values():
            if "294K" in reaction.xs:
                self.assertIn("600K", reaction.xs, f"Missing merged cross section for MT{reaction.mt}")
                np.testing.assert_array_equal(reaction.xs["294K"].y, reaction.xs["600K"].y)
        self.assert_catalog(xml, {"He3"})

    def test_duplicate_rounded_temperatures_fail_before_output(self):
        second = self.work / "duplicate.ace"
        second.write_bytes(self.he3.read_bytes())
        target = self.work / "out"
        with self.assertRaisesRegex(ValueError, "Duplicate rounded temperature"):
            convert_ace_to_openmc_hdf5([self.he3, second], target)
        self.assertFalse(target.exists())

    def test_output_collisions_require_overwrite_and_catalog_is_current(self):
        target = self.work / "out"
        paths, xml = convert_ace_to_openmc_hdf5(self.he3, target)
        before = {p: p.read_bytes() for p in [*paths, xml]}
        with self.assertRaises(FileExistsError):
            convert_ace_to_openmc_hdf5(self.he3, target)
        for p, content in before.items():
            self.assertEqual(p.read_bytes(), content)
        paths, xml = convert_ace_to_openmc_hdf5(self.u235, target, overwrite=True)
        self.assertEqual([p.name for p in paths], ["U235.h5"])
        self.assert_catalog(xml, {"U235"})

    def test_single_nuclide_custom_h5_path(self):
        target = self.work / "custom" / "chosen.h5"
        paths, xml = convert_ace_to_openmc_hdf5(self.he3, self.work / "unused", h5_path=target)
        self.assertEqual(paths, [target.resolve()])
        self.assertEqual(xml, target.parent / "cross_sections.xml")
        self.assert_catalog(xml, {"He3"})

    def test_custom_h5_path_rejects_multiple_nuclides(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            convert_ace_to_openmc_hdf5([self.he3, self.u235], self.work, h5_path=self.work / "both.h5")
        self.assertFalse((self.work / "both.h5").exists())

    def test_staged_export_failure_does_not_publish_partial_results(self):
        original = self.data_api.IncidentNeutron.export_to_hdf5

        def fail_for_u235(neutron, path, **kwargs):
            if neutron.name == "U235":
                raise OSError("simulated export failure")
            return original(neutron, path, **kwargs)

        target = self.work / "out"
        with patch.object(self.data_api.IncidentNeutron, "export_to_hdf5", fail_for_u235):
            with self.assertRaisesRegex(OSError, "simulated export"):
                convert_ace_to_openmc_hdf5([self.he3, self.u235], target)
        self.assertEqual(list(target.iterdir()), [])

    def test_multigroup_table_after_neutron_table_is_not_silently_skipped(self):
        combined = self.work / "mixed.ace"
        mg = ROOT / "examples" / "mg08" / "99902.02m"
        combined.write_bytes(self.he3.read_bytes().rstrip() + b"\n" + mg.read_bytes())
        target = self.work / "out"
        with self.assertRaises(ValueError):
            convert_ace_to_openmc_hdf5(combined, target)
        self.assertFalse(target.exists())

    def test_truncated_binary_cli_reports_error_without_traceback(self):
        binary = self.work / "truncated.bin"
        binary.write_bytes(b"1001.80c  " + b"\x00" * 39 + b"\xff")
        result = subprocess.run(
            [sys.executable, "-m", "multi_group_ace.cli", "--ace", str(binary), "--outdir", str(self.work / "out")],
            cwd=ROOT, capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_invalid_temperature_metadata_is_rejected_before_output(self):
        for temperature in ("-2.5301E-08", "nan", "inf"):
            with self.subTest(temperature=temperature):
                lines = self.he3.read_text().splitlines(keepends=True)
                lines[1] = lines[1].replace("2.5301E-08", temperature)
                source = self.work / "bad-temperature.ace"
                source.write_text("".join(lines))
                target = self.work / "out"
                with self.assertRaises(ValueError):
                    convert_ace_to_openmc_hdf5(source, target)
                self.assertFalse(target.exists())

    def test_invalid_atomic_weight_metadata_is_rejected_before_output(self):
        for ratio in ("0", "-1", "nan", "inf"):
            with self.subTest(ratio=ratio):
                lines = self.he3.read_text().splitlines(keepends=True)
                lines[1] = lines[1].replace("2.989032", ratio)
                source = self.work / "bad-atomic-weight.ace"
                source.write_text("".join(lines))
                target = self.work / "out"
                with self.assertRaises(ValueError):
                    convert_ace_to_openmc_hdf5(source, target)
                self.assertFalse(target.exists())


if __name__ == "__main__":
    unittest.main()
