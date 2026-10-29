# Multi-Group-ACE

Python utilities for multigroup workflows and continuous-energy neutron data conversion:

- build simple multigroup ACE tables from JSON inputs
- convert multigroup ACE tables to OpenMC-style `mgxs.h5` libraries
- convert continuous-energy neutron ACE libraries to OpenMC nuclide HDF5 files and `cross_sections.xml`
- generate `geom.msh` meshes from `.geo` / `.step` files or supported MCNP surface cards
- prepare OpenSn-facing assets from bundled verification-style problem directories

This repository was developed as an undergraduate research project aimed at making it easier to move from MCNP-style multigroup data and geometry inputs toward assets that can be used in OpenSn workflows. It is centered on the `multi_group_ace` package and a set of example cases under `examples/`.

The `--ace` option detects the table type from file contents. Legacy ASCII multigroup tables (names ending in `m`) produce `mgxs.h5`; continuous-energy neutron tables (names ending in `c`, including `nc`) produce OpenMC continuous-energy HDF5. The bundled `examples/` tables are multigroup data and cannot recover continuous-energy detail. See [continuous-energy conversion](#continuous-energy-ace-for-openmc) for the new workflow.

## Project Goal

The long-term goal of this repository is to help simplify the path from MCNP-related problem data to OpenSn-ready inputs and supporting files. The current code already covers several important pieces of that workflow:

- defining multigroup data in JSON and writing ACE tables
- reading ACE tables and converting them into `mgxs.h5`
- converting existing continuous-energy neutron ACE evaluations into OpenMC's native data format
- generating Gmsh meshes from explicit geometry files or supported MCNP surface cards
- organizing verification-style problem outputs in a way that can support future OpenSn use

## Current Status

This repository is best viewed as a working research prototype rather than a finished end-to-end production tool.

- The package generates multigroup ACE tables, `mgxs.h5`, and `geom.msh`, and converts continuous-energy neutron ACE data using OpenMC's Python API.
- The repository includes a Python OpenSn driver example in `examples/mg08_out/run_keigen.py`.
- The project does **not** yet generate a complete standalone OpenSn text input deck automatically.
- Some workflows are currently easier to describe as "OpenSn-facing" than "fully automated OpenSn execution."

## Quick Start

Run commands from the repository root.

Install the optional dependencies used by the main workflows:

```bash
python -m pip install h5py gmsh
```

Generate an ACE table from a JSON material description:

```bash
python -m multi_group_ace.cli --config "json examples/ex4.json"
```

Convert a multigroup ACE table to `mgxs.h5`:

```bash
python -m multi_group_ace.cli --ace "examples/mg08/99902.02m" --h5 "examples/mg08_out/mgxs.h5" --library-name 99902.02m
```

For continuous-energy conversion, install the optional OpenMC API and pass one or more neutron ACE files:

```bash
python -m pip install -r requirements-ce.txt
python -m multi_group_ace.cli --ace path/to/He3.ace path/to/U235.ace --outdir openmc_data
```

This writes `openmc_data/He3.h5`, `openmc_data/U235.h5`, and `openmc_data/cross_sections.xml`. The filenames shown assume those nuclides are present in the input tables. An OpenMC transport executable is not required for conversion.

Batch-convert a verification tree and generate meshes:

```bash
python -m multi_group_ace.cli --verification-root examples --outdir build/verification --make-mesh --gmsh gmsh
```

## What the Repository Produces

Depending on the workflow, the package writes some combination of:

- simple multigroup ACE tables
- `xsdir_mcnp6.3`
- summary text files such as `<base>.summary.txt`
- OpenSn-compatible multigroup HDF5 libraries (`mgxs.h5`)
- OpenMC continuous-energy neutron HDF5 files (`<nuclide>.h5`) and an XML cross-section catalog
- Gmsh meshes (`geom.msh`)

For OpenSn-related usage, the repository currently prepares solver-facing inputs such as `mgxs.h5` and `geom.msh`, but it does **not** yet generate a standalone OpenSn text input deck.

## Workflow at a Glance

```text
JSON material description
    -> MGData
    -> ACE arrays (NXS / JXS / XSS)
    -> ACE file

Multigroup ACE file
    -> MGData
    -> mgxs.h5

Continuous-energy neutron ACE file(s)
    -> OpenMC IncidentNeutron data, merged by nuclide and temperature
    -> <nuclide>.h5 + cross_sections.xml

.geo / .step / .stp / MCNP input
    -> Gmsh geometry handling
    -> geom.msh

mgxs.h5 + geom.msh
    -> examples/mg08_out/run_keigen.py
    -> OpenSn solve
```

## Repository Layout

| Path | Role |
| --- | --- |
| `multi_group_ace/` | Main package implementation |
| `examples/` | Verification-style source cases and generated outputs |
| `json examples/` | Sample JSON material definitions |
| `doc/` | LaTeX and PDF documentation |
| `notes/` | Project notes and reference text |

Core modules inside `multi_group_ace/`:

- `core.py`: shared constants and the `MGData` container
- `json_loader.py`: JSON parsing, validation, and `MGData` construction
- `converter.py`: `MGData` to ACE `NXS`, `JXS`, `XSS` construction
- `ace_io.py`: ACE reading/writing plus summary and `xsdir` output
- `hdf5_io.py`: OpenMC-style multigroup HDF5 writing
- `continuous_energy.py`: content detection, OpenMC neutron ACE conversion, temperature merging, and XML catalog generation
- `gmsh_geometry.py`: geometry handling and MCNP-to-Gmsh mesh generation
- `cli.py`: command-line entry point

## Main Interface

The primary interface is:

```bash
python -m multi_group_ace.cli --help
```

Available arguments:

- `--config CONFIG`: JSON input for JSON-to-ACE conversion
- `--ace FILE [FILE ...]`: one multigroup ACE input, or one or more continuous-energy neutron ACE libraries; routing uses table contents rather than filename extensions
- `--h5 H5`: output HDF5 path; a continuous-energy override must end in `.h5` and produce exactly one nuclide
- `--library-name LIBRARY_NAME`: optional dataset/group name inside `mgxs.h5` (single-file `--ace` mode only)
- `--verification-root VERIFICATION_ROOT`: root directory for recursive multigroup ACE scanning
- `--outdir OUTDIR`: output root for scanned multigroup libraries or continuous-energy conversion (default: `h5_out`)
- `--metastable-scheme {mcnp,nndc}`: continuous-energy ZAID interpretation (default: `mcnp`)
- `--overwrite`: replace existing continuous-energy HDF5 files and XML catalog
- `--make-mesh`: also generate geometry meshes during recursive scan mode
- `--gmsh GMSH`: path to the Gmsh executable

Choose exactly one of `--config`, `--ace`, and `--verification-root`. The historical script entry point, `python simple_ace_mg_full.py`, delegates to the same maintained package CLI; its old helper implementations are retained for reference.

`--library-name` is multigroup-only. `--overwrite` and `--metastable-scheme` are continuous-energy-only. A continuous-energy `--h5` override cannot be combined with an explicit `--outdir`. Recursive verification scanning remains a multigroup workflow.

## Core Workflows

### 1. JSON to ACE

Use a JSON material definition to build a simple multigroup ACE table.

```bash
python -m multi_group_ace.cli --config "json examples/ex4.json"
```

Typical outputs:

- ACE table file
- `xsdir_mcnp6.3`
- `<base>.summary.txt`

Expected JSON fields include:

- `zaid`
- `NG`
- `erg`
- `sigt`
- `sigc`
- `sigf`
- `nu`
- `chi`
- `sigs`

Optional fields include `sig1`, `is_micro`, `atom_density`, `comment`, `tmp`, `nbins`, and `file`.

### 2. ACE to `mgxs.h5`

Convert a multigroup ACE table into an OpenSn-compatible multigroup HDF5 library.

```bash
python -m multi_group_ace.cli --ace "examples/mg08/99902.02m" --h5 "examples/mg08_out/mgxs.h5" --library-name 99902.02m
```

The resulting file contains:

- root-level MGXS metadata
- energy-group structure stored in ascending eV
- one or more material groups
- total, absorption, fission, chi, nu-fission, and scattering datasets

### 3. Recursive Verification Conversion

Scan a directory tree, identify multigroup ACE tables, group them by problem directory, and write one combined `mgxs.h5` per problem.

```bash
python -m multi_group_ace.cli --verification-root examples --outdir build/verification
```

To also generate meshes:

```bash
python -m multi_group_ace.cli --verification-root examples --outdir build/verification --make-mesh --gmsh gmsh
```

## Geometry and Meshing

When mesh generation is enabled, the code follows this order:

1. If a problem directory contains a `.geo`, `.step`, or `.stp` file, the Gmsh executable is called directly on that geometry.
2. Otherwise, the package searches for an MCNP input file ending in `.inp`, `.i`, or `.mcnp` and parses supported surface cards.

Supported MCNP surface types include:

- planes: `px`, `py`, `pz`
- spheres: `so`, `s`, `sph`
- cylinders: `cz`, `cx`, `cy`, `rcc`

The generated mesh is written as `geom.msh`. Boundary surfaces may also be tagged with names such as `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, and `zmax` to make downstream boundary-condition setup easier.

The MCNP fallback builds simplified geometry from surface cards. It does not translate general cell Boolean expressions, material assignments, transformations, or reflecting/periodic boundary conditions. Slab transverse dimensions and infinite-domain proxies need review before a solver run, and sphere/cylinder builders assume concentric regions. A mesh being readable does not establish that it reproduces the original transport problem. Prefer an explicit `.geo`/STEP geometry when exact geometry is required.

## Multigroup HDF5 Output Format

The HDF5 writer creates an OpenSn-compatible MGXS structure. At a high level:

```text
/
  attrs: filetype, version, energy_groups, delayed_groups, group structure
  /<material_name>
    attrs: atomic_weight_ratio, fissionable, representation, scatter_format, scatter_shape, order
    /kTs
    /<temperature_group>
      /total
      /absorption
      /fission
      /chi
      /nu-fission
      /scatter_data
```

If anisotropic scattering is present, the file is written with scattering order 1 and stores both P0 and P1 terms in the flattened scatter matrix. Otherwise it is written as isotropic with scattering order 0.

The `representation` attribute remains `isotropic` for both P0 and P1: it describes dependence on the incoming direction, while `order` describes scattering anisotropy. P1 recovered from ACE equiprobable angular bins is approximate; higher angular moments are not preserved. Combined libraries require matching energy boundaries and unique table names, and each material retains its own temperature.

## Continuous-Energy ACE for OpenMC

This workflow delegates ACE parsing and neutron HDF5 serialization to OpenMC's [IncidentNeutron API](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.data.IncidentNeutron.html). It processes every table in each supplied ACE library, accepts the ASCII and binary neutron formats supported by OpenMC, and groups compatible temperatures into one GNDS-named file per nuclide. OpenMC is imported only when this workflow is requested.

Install the pinned OpenMC 0.16.0 Python API, then convert files from the same evaluation:

```bash
python -m pip install -r requirements-ce.txt
python -m multi_group_ace.cli --ace U235_294K.ace U235_600K.ace --outdir openmc_data --metastable-scheme mcnp
```

That example produces `U235.h5` containing both temperatures and a `cross_sections.xml` catalog. For a single resulting nuclide, a custom output filename is also supported:

```bash
python -m multi_group_ace.cli --ace U235.ace --h5 openmc_data/U235.h5
```

The custom filename must end in `.h5`; its parent receives `cross_sections.xml`. By default, any existing destination HDF5 file or catalog causes an error before export. Use `--overwrite` to replace them. The catalog describes only the current invocation's outputs; overwrite does not append entries from an older catalog.

Use the ZAID metastable convention matching the source library (`mcnp`, the default, or `nndc`). Inputs merged for a nuclide must come from the same evaluation and have compatible reaction sets and atomic weight ratios. Duplicate rounded temperature labels and incompatible inputs are rejected rather than silently dropping data. Thermal scattering tables, other particle types, and mixed multigroup/continuous-energy inputs are unsupported.

The public API returns absolute paths for the HDF5 outputs and catalog:

```python
from multi_group_ace import convert_ace_to_openmc_hdf5

hdf5_paths, xml_path = convert_ace_to_openmc_hdf5(
    ["He3.ace", "U235.ace"], "openmc_data",
    metastable_scheme="mcnp", overwrite=False,
)
```

Point OpenMC at the generated catalog, for example in PowerShell:

```powershell
$env:OPENMC_CROSS_SECTIONS = (Resolve-Path "openmc_data/cross_sections.xml").Path
```

The neutron files use `filetype=data_neutron` with energy grids and reaction data; they are different from `filetype=mgxs` libraries. Conversion does not collapse continuous-energy data into groups or reconstruct continuous-energy information from multigroup inputs. Nuclides and reactions needed for a transport model must all be present in its data library.

## Verification

Run the regression suite from the repository root:

```bash
python -m unittest discover -s tests -v
```

The original 24 regression checks cover JSON/ACE/MGXS workflows, all 15 bundled multigroup ACE tables across 13 problems, neutron balance, five infinite-medium analytical benchmarks, angular bins, invalid inputs, and generation/reopening of all 13 example meshes. Additional tests cover continuous-energy routing, validation, and real OpenMC round trips. Tests write generated output into temporary directories.

For the complete continuous-energy checks on Windows, use an isolated environment and fetch the checksum-verified neutron ACE fixtures:

```powershell
python -m venv --system-site-packages .venv-ce
.venv-ce\Scripts\python -m pip install -r requirements-ce.txt
.venv-ce\Scripts\python tests/fetch_ce_data.py
$env:MULTI_GROUP_ACE_REQUIRE_CE_TESTS = "1"
.venv-ce\Scripts\python -m unittest discover -s tests -v
```

The fixture script downloads authentic He3 and U235 ENDF/B-VII.1 tables from the NJOY ACEtk test resources into ignored `.test-data/ce/` storage and records their provenance. Real integration tests compare conversion outputs with direct OpenMC ACE loading, reopen HDF5 and the XML catalog, and exercise multiple tables, binary inputs, and temperature merging. Modified-temperature copies are labeled synthetic plumbing tests. Optional checks skip when dependencies or fixtures are absent; a release verification run must have no skips.

The complete verification run passed **51 tests with no skips or failures**: 24 original workflow checks, 14 continuous-energy unit checks, and 13 real OpenMC integration checks. The required-CE environment setting above makes missing OpenMC or fixture prerequisites an error.

These checks use OpenMC's Python data API and do not execute its transport solver or OpenSn, or establish finite-geometry benchmark agreement. See [the repository audit](doc/repository_audit.md) for measured results and remaining limits.

## Examples

The `examples/` directory contains bundled source cases and generated outputs. Typical contents include:

- `mgXX/description.json`
- `mgXX/*.inp`
- `mgXX/*ACE table*`
- `mgXX_out/mgxs.h5`
- `mgXX_out/geom.msh`

The most OpenSn-facing example in the repository is:

```bash
python examples/mg08_out/run_keigen.py
```

That script loads `geom.msh` and `mgxs.h5` through `pyopensn`, constructs a discrete ordinates problem, runs a k-eigenvalue solve, and writes ParaView-compatible output.

For presentation purposes, this example is best understood as a proof of concept showing how the repository's generated assets can be used downstream in OpenSn.

## Requirements and Optional Dependencies

The base package uses the Python standard library plus optional dependencies depending on the workflow:

- `h5py`: required for writing `mgxs.h5`
- `openmc`: optional for continuous-energy conversion; `requirements-ce.txt` pins the tested Python API to 0.16.0 and installs its Python dependencies
- `gmsh` Python package: required when converting MCNP input text directly through the Gmsh Python API
- Gmsh executable: used when meshing `.geo`, `.step`, or `.stp` geometry files
- `pyopensn`: only needed for `examples/mg08_out/run_keigen.py`

Example installation:

```bash
python -m pip install h5py gmsh
```

If you want to run the OpenSn example driver, install your local `pyopensn` environment separately.

The pinned OpenMC API supports Python 3.12 or newer and installs from its official Git source tag, so Git and network access are needed for that dependency installation. The multigroup package itself does not import OpenMC at startup.

## Notes

- Internally, scattering data are stored with `row = destination group` and `column = source group`.
- External formats such as ACE and HDF5 are reordered as needed when written.
- Batch HDF5 conversion writes one combined `mgxs.h5` per problem directory, not one per ACE file.
- Continuous-energy conversion writes one HDF5 file per nuclide, combining compatible temperatures, plus `cross_sections.xml`.
- The repository currently focuses on preparing data and mesh assets; the final OpenSn automation layer is still a future direction for the project.

## Documentation

Additional documentation is available in `doc/`:

- [`doc/multi_group_ace_technical_documentation.tex`](doc/multi_group_ace_technical_documentation.tex)
- [`doc/Multi_Group_ACE_Documentation.pdf`](doc/Multi_Group_ACE_Documentation.pdf)
- [`doc/repository_audit.md`](doc/repository_audit.md)
