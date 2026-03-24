# Multi-Group-ACE

Python utilities for multigroup nuclear data workflows:

- build simple multigroup ACE tables from JSON inputs
- convert ACE tables to OpenMC-style `mgxs.h5` libraries
- generate `geom.msh` meshes from `.geo` / `.step` files or supported MCNP surface cards
- prepare OpenSn-facing assets from bundled verification-style problem directories

This repository was developed as an undergraduate research project aimed at making it easier to move from MCNP-style multigroup data and geometry inputs toward assets that can be used in OpenSn workflows. It is centered on the `multi_group_ace` package and a set of example cases under `examples/`.

## Project Goal

The long-term goal of this repository is to help simplify the path from MCNP-related problem data to OpenSn-ready inputs and supporting files. The current code already covers several important pieces of that workflow:

- defining multigroup data in JSON and writing ACE tables
- reading ACE tables and converting them into `mgxs.h5`
- generating Gmsh meshes from explicit geometry files or supported MCNP surface cards
- organizing verification-style problem outputs in a way that can support future OpenSn use

## Current Status

This repository is best viewed as a working research prototype rather than a finished end-to-end production tool.

- The package can already generate useful intermediate assets such as ACE tables, `mgxs.h5`, and `geom.msh`.
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

Convert an ACE table to `mgxs.h5`:

```bash
python -m multi_group_ace.cli --ace "examples/mg08/99902.02m" --h5 "examples/mg08_out/mgxs.h5" --library-name 99902.02m
```

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
- Gmsh meshes (`geom.msh`)

For OpenSn-related usage, the repository currently prepares solver-facing inputs such as `mgxs.h5` and `geom.msh`, but it does **not** yet generate a standalone OpenSn text input deck.

## Workflow at a Glance

```text
JSON material description
    -> MGData
    -> ACE arrays (NXS / JXS / XSS)
    -> ACE file

ACE file
    -> MGData
    -> mgxs.h5

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
- `gmsh_geometry.py`: geometry handling and MCNP-to-Gmsh mesh generation
- `cli.py`: command-line entry point

## Main Interface

The primary interface is:

```bash
python -m multi_group_ace.cli --help
```

Available arguments:

- `--config CONFIG`: JSON input for JSON-to-ACE conversion
- `--ace ACE`: ACE input for ACE-to-HDF5 conversion
- `--h5 H5`: output HDF5 path
- `--library-name LIBRARY_NAME`: optional dataset/group name inside `mgxs.h5`
- `--verification-root VERIFICATION_ROOT`: root directory for recursive ACE scanning
- `--outdir OUTDIR`: output root for scanned HDF5 files
- `--make-mesh`: also generate geometry meshes during recursive scan mode
- `--gmsh GMSH`: path to the Gmsh executable

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

Convert an ACE table into an OpenSn-compatible multigroup HDF5 library.

```bash
python -m multi_group_ace.cli --ace "examples/mg08/99902.02m" --h5 "examples/mg08_out/mgxs.h5" --library-name 99902.02m
```

The resulting file contains:

- root-level MGXS metadata
- energy-group structure stored in ascending eV
- one or more material groups
- total, absorption, fission, chi, nu-fission, and scattering datasets

### 3. Recursive Verification Conversion

Scan a directory tree, identify ACE tables, group them by problem directory, and write one combined `mgxs.h5` per problem.

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
2. Otherwise, the package searches for an MCNP input file such as `.inp`, `.i`, or `.txt` and parses supported surface cards.

Supported MCNP surface types include:

- planes: `px`, `py`, `pz`
- spheres: `so`, `s`, `sph`
- cylinders: `cz`, `cx`, `cy`, `rcc`

The generated mesh is written as `geom.msh`. Boundary surfaces may also be tagged with names such as `xmin`, `xmax`, `ymin`, `ymax`, `zmin`, and `zmax` to make downstream boundary-condition setup easier.

## HDF5 Output Format

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
- `gmsh` Python package: required when converting MCNP input text directly through the Gmsh Python API
- Gmsh executable: used when meshing `.geo`, `.step`, or `.stp` geometry files
- `pyopensn`: only needed for `examples/mg08_out/run_keigen.py`

Example installation:

```bash
python -m pip install h5py gmsh
```

If you want to run the OpenSn example driver, install your local `pyopensn` environment separately.

## Notes

- Internally, scattering data are stored with `row = destination group` and `column = source group`.
- External formats such as ACE and HDF5 are reordered as needed when written.
- Batch HDF5 conversion writes one combined `mgxs.h5` per problem directory, not one per ACE file.
- The repository currently focuses on preparing data and mesh assets; the final OpenSn automation layer is still a future direction for the project.

## Documentation

Additional documentation is available in `doc/`:

- [`doc/multi_group_ace_technical_documentation.tex`](doc/multi_group_ace_technical_documentation.tex)
