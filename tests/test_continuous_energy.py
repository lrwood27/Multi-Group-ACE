"""Dependency-independent routing, validation, and CLI regression checks."""

import contextlib
import io
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from multi_group_ace import continuous_energy as ce
from multi_group_ace.ace_io import read_mg_ace
from multi_group_ace.cli import main


ROOT = Path(__file__).resolve().parents[1]
MG = ROOT / "examples" / "mg08" / "99902.02m"


class ContinuousEnergyChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def fixture(self, content, name="input.data"):
        path = self.work / name
        path.write_bytes(content)
        return path

    def cli_error(self, *arguments):
        output = io.StringIO()
        with contextlib.redirect_stderr(output), self.assertRaises(SystemExit) as result:
            main(list(map(str, arguments)))
        self.assertEqual(result.exception.code, 2)
        self.assertNotIn("Traceback", output.getvalue())
        return output.getvalue()

    def test_detection_uses_contents_instead_of_filename(self):
        neutron = self.fixture(b"1001.80c 1.0 2.5301e-8\n", "misleading.01m")
        self.assertEqual(ce.detect_ace_type(neutron), "continuous-energy")
        multigroup = self.fixture(MG.read_bytes(), "misleading.80c")
        self.assertEqual(ce.detect_ace_type(multigroup), "multigroup")

    def test_detection_handles_new_ascii_and_multiletter_table_suffixes(self):
        for version, table in (("2.0.0", "2003.710nc"), ("2.0.1", "92235.710nc")):
            with self.subTest(version=version):
                path = self.fixture(f"{version} {table} ENDFB-VII.1\n".encode())
                self.assertEqual(ce.detect_ace_type(path), "continuous-energy")

    def test_detection_handles_binary_and_unsupported_tables(self):
        binary = self.fixture(b"1001.80c  " + b"\x00" * 100)
        self.assertEqual(ce.detect_ace_type(binary), "continuous-energy")
        for name in ("lwtr.20t", "1001.80p"):
            with self.subTest(name=name):
                path = self.fixture(f"{name} 1.0 2.5301e-8\n".encode())
                self.assertEqual(ce.detect_ace_type(path), "unsupported")

    def test_empty_and_unrecognized_headers_fail(self):
        for content in (b"", b"not ACE data\n"):
            with self.subTest(content=content), self.assertRaises(ValueError):
                ce.detect_ace_type(self.fixture(content))

    def test_missing_openmc_has_install_guidance(self):
        with patch.dict(sys.modules, {"openmc": None, "openmc.data": None}):
            with self.assertRaisesRegex(ImportError, "requirements-ce.txt"):
                ce._load_openmc()

    def test_multigroup_and_invalid_input_are_rejected_before_openmc_import(self):
        unsupported = self.fixture(b"lwtr.20t 1.0 2.5301e-8\n")
        with patch.object(ce, "_load_openmc", side_effect=AssertionError("should not import")):
            for paths in ([], [MG], [unsupported], [MG, MG]):
                with self.subTest(paths=paths), self.assertRaises(ValueError):
                    ce.convert_ace_to_openmc_hdf5(paths, self.work / "out")

    def test_invalid_metastable_scheme_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "metastable_scheme"):
            ce.convert_ace_to_openmc_hdf5([MG], self.work, metastable_scheme="invalid")

    def test_incompatible_reaction_sets_are_rejected_before_export(self):
        source = self.fixture(b"1001.80c 1.0 2.5301e-8\n")
        tables = [
            SimpleNamespace(name="1001.80c", temperature=2.53e-8, atomic_weight_ratio=1.),
            SimpleNamespace(name="1001.81c", temperature=5.17e-8, atomic_weight_ratio=1.),
        ]
        models = [
            SimpleNamespace(name="H1", kTs=[0.0253], temperatures=["294K"], atomic_weight_ratio=1., reactions={2: object(), 102: object()}),
            SimpleNamespace(name="H1", kTs=[0.0517], temperatures=["600K"], atomic_weight_ratio=1., reactions={2: object()}),
        ]
        api = SimpleNamespace(
            ace=SimpleNamespace(Library=lambda path: SimpleNamespace(tables=tables)),
            IncidentNeutron=SimpleNamespace(from_ace=lambda table, **kwargs: models[tables.index(table)]),
        )
        with patch.object(ce, "_load_openmc", return_value=api):
            with self.assertRaisesRegex(ValueError, "Incompatible reaction sets"):
                ce.convert_ace_to_openmc_hdf5(source, self.work / "out")
        self.assertFalse((self.work / "out").exists())

    def test_cli_rejects_mixed_input_files(self):
        neutron = self.fixture(b"1001.80c 1.0 2.5301e-8\n")
        self.assertIn("Cannot mix", self.cli_error("--ace", MG, neutron))

    def test_cli_rejects_multiple_multigroup_inputs(self):
        self.assertIn("one input", self.cli_error("--ace", MG, MG))

    def test_cli_rejects_continuous_energy_only_options_on_multigroup_input(self):
        for options in (("--overwrite",), ("--metastable-scheme", "nndc")):
            with self.subTest(options=options):
                self.assertIn("continuous-energy-only", self.cli_error("--ace", MG, *options))

    def test_cli_rejects_ambiguous_or_multigroup_only_ce_options(self):
        neutron = self.fixture(b"1001.80c 1.0 2.5301e-8\n")
        self.assertIn("cannot be combined", self.cli_error("--ace", neutron, "--h5", "out.h5", "--outdir", "out"))
        self.assertIn("multigroup-only", self.cli_error("--ace", neutron, "--library-name", "custom"))

    def test_cli_rejects_thermal_scattering(self):
        thermal = self.fixture(b"lwtr.20t 1.0 2.5301e-8\n")
        self.assertIn("Thermal scattering", self.cli_error("--ace", thermal))

    def test_multigroup_first_concatenated_input_does_not_ignore_trailing_table(self):
        mixed = self.fixture(MG.read_bytes().rstrip() + b"\n1001.80c 1.0 2.5301e-8\n")
        with self.assertRaises(ValueError):
            read_mg_ace(mixed)
        self.cli_error("--ace", mixed, "--h5", self.work / "out.h5")
        self.assertFalse((self.work / "out.h5").exists())


if __name__ == "__main__":
    unittest.main()
