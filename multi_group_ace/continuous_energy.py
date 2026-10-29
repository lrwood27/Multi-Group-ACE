"""Convert continuous-energy neutron ACE libraries using OpenMC's data API.

OpenMC is optional and is imported only when a continuous-energy conversion
is requested. Multigroup ACE is handled separately by :mod:`ace_io`.
"""

from collections import defaultdict
from copy import deepcopy
from math import isclose, isfinite
import os
from pathlib import Path
import re
from struct import error as StructError
from tempfile import TemporaryDirectory


_TABLE_NAME = re.compile(r"[A-Za-z0-9_-]+\.\d+[A-Za-z]*([A-Za-z])\Z")


def detect_ace_type(path):
    """Identify the first ACE table from its contents, without importing OpenMC.

    Return ``multigroup``, ``continuous-energy``, or ``unsupported``. The
    filesystem extension is irrelevant. Both ASCII header generations and
    the fixed-width name at the start of a binary ACE header are recognized.
    This is routing only; the relevant reader validates the entire table.
    """
    with Path(path).open("rb") as stream:
        header = stream.read(512)
    if not header.strip():
        raise ValueError(f"Empty ACE input: {path}")

    candidates = []
    words = header.splitlines()[0].decode("ascii", errors="replace").split()
    if words:
        candidates.append(words[0])
        if re.fullmatch(r"\d+\.\d+\.\d+", words[0]) and len(words) > 1:
            candidates.append(words[1])
    candidates.append(header[:10].decode("ascii", errors="replace").strip())
    for name in candidates:
        match = _TABLE_NAME.fullmatch(name)
        if match:
            return {"m": "multigroup", "c": "continuous-energy"}.get(
                match[1].lower(), "unsupported"
            )
    raise ValueError(f"Unrecognized ACE table header in {path}.")


def _load_openmc():
    try:
        import openmc.data
    except ImportError as exc:
        raise ImportError(
            "Continuous-energy conversion requires the OpenMC Python API. "
            "Install it with: python -m pip install -r requirements-ce.txt "
            f"(import failed: {exc})"
        ) from exc
    return openmc.data


def convert_ace_to_openmc_hdf5(
    ace_files, output_dir, *, h5_path=None, metastable_scheme="mcnp", overwrite=False
):
    """Write neutron HDF5 files and a ``cross_sections.xml`` catalog.

    ``ace_files`` is a path or iterable of paths. Every table must contain
    continuous-energy neutron data. Compatible temperatures from the same
    evaluation are merged into one GNDS-named file per nuclide. Duplicate
    rounded temperature labels and mismatched reaction sets are rejected.

    ``h5_path`` overrides the output location for a single resulting nuclide;
    its parent also receives the XML catalog. Otherwise ``output_dir`` gets
    all outputs. Existing destinations require ``overwrite=True``, which
    replaces the catalog with entries for this invocation only.

    Return ``(hdf5_paths, xml_path)`` as absolute :class:`pathlib.Path` objects.
    OpenMC's Python API is required; its transport executable is not.
    """
    if metastable_scheme not in ("mcnp", "nndc"):
        raise ValueError("metastable_scheme must be 'mcnp' or 'nndc'.")
    if isinstance(ace_files, (str, os.PathLike)):
        ace_files = [ace_files]
    sources = [Path(path).resolve() for path in ace_files]
    if not sources:
        raise ValueError("Provide at least one continuous-energy neutron ACE file.")
    if len(sources) != len(set(sources)):
        raise ValueError("Duplicate ACE input paths are not allowed.")
    for source in sources:
        kind = detect_ace_type(source)
        if kind == "multigroup":
            raise ValueError(
                f"{source} contains multigroup ACE data. Continuous-energy "
                "detail cannot be recovered from multigroup cross sections."
            )
        if kind != "continuous-energy":
            raise ValueError(
                f"{source} is unsupported: only continuous-energy neutron ACE "
                "tables (names ending in 'c') can be converted; thermal "
                "scattering and other table types are not supported."
            )

    data_api = _load_openmc()
    grouped = defaultdict(list)
    for source in sources:
        try:
            tables = data_api.ace.Library(str(source)).tables
        except (ValueError, IndexError, KeyError, TypeError, AssertionError, StructError) as exc:
            raise ValueError(f"Could not read ACE file {source}: {exc}") from exc
        if not tables:
            raise ValueError(f"No ACE tables found in {source}.")
        for table in tables:
            if not table.name.lower().endswith("c"):
                raise ValueError(
                    f"Unsupported table {table.name} in {source}: every table "
                    "must be continuous-energy neutron data (ending in 'c')."
                )
            if not isfinite(table.temperature) or table.temperature < 0:
                raise ValueError(f"Invalid temperature in table {table.name}: expected finite, nonnegative kT.")
            if not isfinite(table.atomic_weight_ratio) or table.atomic_weight_ratio <= 0:
                raise ValueError(f"Invalid atomic weight ratio in table {table.name}: expected a finite positive value.")
            try:
                # OpenMC converts some XSS slices in place (notably heating
                # from MeV to eV). Preserve pristine tables for temperature
                # merging, which parses them again through the same API.
                neutron = data_api.IncidentNeutron.from_ace(
                    deepcopy(table), metastable_scheme=metastable_scheme
                )
            except (ValueError, IndexError, KeyError, TypeError, AssertionError) as exc:
                raise ValueError(f"Could not convert table {table.name} in {source}: {exc}") from exc
            grouped[neutron.name].append((table, neutron, source))

    # OpenMC otherwise warns and skips duplicate temperatures/new reactions.
    # Validate all inputs before merging so data cannot disappear silently.
    models = {}
    for name in sorted(grouped):
        entries = sorted(grouped[name], key=lambda entry: entry[1].kTs[0])
        reference = entries[0][1]
        seen_temperatures = set()
        for table, neutron, source in entries:
            temperature = neutron.temperatures[0]
            if temperature in seen_temperatures:
                raise ValueError(
                    f"Duplicate rounded temperature {temperature} for {name} "
                    f"(table {table.name} in {source})."
                )
            seen_temperatures.add(temperature)
            if set(neutron.reactions) != set(reference.reactions):
                raise ValueError(f"Incompatible reaction sets for {name}; merge temperatures from the same evaluation.")
            if not isclose(neutron.atomic_weight_ratio, reference.atomic_weight_ratio, rel_tol=1e-8):
                raise ValueError(f"Incompatible atomic weight ratios for {name}.")
        for table, _, _ in entries[1:]:
            reference.add_temperature_from_ace(table, metastable_scheme=metastable_scheme)
        models[name] = reference

    if h5_path is not None:
        if len(models) != 1:
            raise ValueError("--h5/h5_path requires exactly one resulting nuclide; use --outdir instead.")
        destination = Path(h5_path).resolve()
        if destination.suffix != ".h5":
            raise ValueError("The continuous-energy HDF5 output path must end in '.h5'.")
        directory = destination.parent
        hdf5_paths = [destination]
    else:
        directory = Path(output_dir).resolve()
        hdf5_paths = [directory / f"{name}.h5" for name in models]
    xml_path = directory / "cross_sections.xml"
    destinations = [*hdf5_paths, xml_path]
    for destination in destinations:
        if destination in sources:
            raise ValueError(f"Output would replace an ACE input: {destination}")
        if destination.exists():
            if not destination.is_file():
                raise ValueError(f"Output destination is not a file: {destination}")
            if not overwrite:
                raise FileExistsError(f"Output already exists: {destination}. Use --overwrite to replace it.")

    directory.mkdir(parents=True, exist_ok=True)
    # Complete conversion and XML generation before publishing any files.
    # Stage on the same filesystem so each final replacement is atomic.
    with TemporaryDirectory(prefix=".ace-conversion-", dir=directory) as staging:
        staging = Path(staging)
        library = data_api.DataLibrary()
        for neutron, destination in zip(models.values(), hdf5_paths):
            staged = staging / destination.name
            neutron.export_to_hdf5(staged, mode="w")
            library.register_file(staged)
            # Registering inspects the staged HDF5, but XML must refer to the
            # permanent filename. DataLibrary has a public list interface.
            library[-1]["path"] = str(destination)
        staged_xml = staging / xml_path.name
        library.export_to_xml(staged_xml)
        for destination in destinations:
            os.replace(staging / destination.name, destination)

    return hdf5_paths, xml_path
