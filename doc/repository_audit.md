# Repository audit — October 29, 2026

The maintained workflows were exercised on Windows with Python 3.13.15, h5py 3.15.1, NumPy 2.3.2, Gmsh 4.14.1, and the OpenMC 0.16.0 Python API in an isolated `.venv-ce` environment. The complete required-CE run passed **51 tests, with zero skips or failures**: 24 original workflow checks, 14 continuous-energy unit checks, and 13 real OpenMC integration checks. Continuous-energy neutron ACE conversion is implemented through OpenMC's own data API; its transport executable and pyopensn are not required for conversion and were not used in this audit.

## Checks performed

- JSON -> ACE -> HDF5 through the CLI, including nested output paths and named material groups.
- Reading, rebuilding, and rereading all 15 bundled ACE tables across 13 problem directories; verification of reaction/scattering balance and data preservation.
- Batch generation of 13 combined HDF5 libraries, checking material names and reaction data.
- Independent zero-leakage eigenvalue calculations for mg01, mg05, mg31, mg74, and mg75; each agrees with the bundled analytical reference to an absolute tolerance of 1e-6. These calculations test the source data and matrix orientation, not an external transport solver.
- P1 round trips, scattering metadata/layout, and inversion of an analytical angular CDF at 10, 100, 1000, and 2000 bins.
- Invalid/truncated inputs, incompatible energy boundaries, duplicate material names, material temperatures, conflicting CLI modes, and exclusion of a nested batch output tree.
- Full batch conversion with mesh generation: 13 HDF5 libraries and 13 meshes. All meshes reopened in Gmsh with nodes, 3D elements, and physical volume groups.
- A separate explicit `.geo` cube was successfully meshed through the installed Windows Gmsh command wrapper.
- Real continuous-energy He3 and U235 ACE inputs were converted into GNDS-named neutron HDF5 files plus `cross_sections.xml`, and reopened with OpenMC's Python API.
- The continuous-energy integration suite checks numerical contents against direct ACE loading, including persisted reaction data, temperatures, energy grids, and available fission/distribution data; it also covers concatenated tables, binary ACE, and synthetic temperature merging.

Run the repeatable regression suite with:

```powershell
.venv-ce\Scripts\python tests/fetch_ce_data.py
$env:MULTI_GROUP_ACE_REQUIRE_CE_TESTS = "1"
.venv-ce\Scripts\python -m unittest discover -s tests -v
```

Install the pinned API first with `.venv-ce\Scripts\python -m pip install -r requirements-ce.txt`. The ordinary Python regression command remains available for multigroup work; integration checks may skip when OpenMC or its downloaded fixtures are absent. Complete release verification requires those checks to run without skips.

## Continuous-energy fixture provenance

The originally planned IAEA FENDL download endpoints returned HTTP 403/Cloudflare responses. Verification therefore uses authentic ENDF/B-VII.1 He3 (`2003.710nc`) and U235 (`92235.710nc`) ACE tables from the primary NJOY ACEtk repository, pinned to commit `aac9e996d5f6c7f5d0729058c0db91cd2f05109e` under `test/resources/`. The fetch script records source URLs and checksums. Downloads stay in ignored `.test-data/ce/` storage and are not bundled into the repository. Temperature-modified copies test parser and merge plumbing only; they are explicitly synthetic and do not represent newly evaluated nuclear data.

## Defects corrected

| Problem | Correction |
| --- | --- |
| Reader could loop forever at EOF in NXS/JXS; batch scanning could consequently hang | Detect unsupported headers and truncated arrays; validate core pointers and finite payload values |
| Angular-bin inversion used the requested bin count as an index into a fixed 1001-point CDF | Search the CDF independently of the requested output count; preserve endpoints and reject negative angular probability densities |
| P1 output claimed incoming-angle-dependent representation without angular dimensions | Keep `representation=isotropic` and record scattering anisotropy with Legendre `order=1` |
| Combined libraries checked energy-grid lengths but not boundaries | Check boundary values before opening the output file |
| All materials inherited the first material's temperature; zero kT was labeled 294K | Compute each material's temperature separately and label zero kT as 0K |
| Duplicate names could fail partway through writing | Reject duplicates before opening the output file |
| Explicit Gmsh geometry used relative paths with a changed working directory; Windows `gmsh.bat` was not resolved | Use absolute input/output paths and resolve the command with `shutil.which` |
| CLI accepted conflicting modes and silently ignored some flags; scanner accepted missing/empty roots | Validate modes/options and scan roots, and exclude the output tree from traversal |
| Historical command-line script ran a separate implementation | Delegate its executable entry point to the maintained package CLI |
| OpenMC's parser changes some ACE XSS slices in place, so parsing a temperature twice could rescale heating by another factor of 1e6 | Use deep copies for preflight conversion and retain pristine tables for temperature merging; real round-trip tests verify preservation |
| Invalid continuous-energy kT or atomic weight metadata could reach OpenMC | Reject negative/nonfinite kT and nonpositive/nonfinite atomic weight ratios before conversion |

## Continuous-energy behavior

`--ace FILE [FILE ...]` detects ACE table types from contents. A single legacy ASCII multigroup table follows the original MGXS path. Continuous-energy neutron libraries use OpenMC's ACE reader to enumerate every table and its `IncidentNeutron` API to serialize neutron HDF5. Compatible temperatures are combined into one file per nuclide; duplicate rounded temperature labels, incompatible reactions or atomic weight ratios, mixed table types, and unsupported thermal/particle tables are rejected.

OpenMC is loaded lazily. Continuous-energy output defaults to `h5_out/`, or the supplied `--outdir`. A `--h5` override must end in `.h5`, must produce one nuclide, and cannot be combined with an explicit `--outdir`. The XML catalog is written beside the outputs using resolvable paths. Existing destination files require `--overwrite`; replacement catalogs contain this invocation's outputs only. The default metastable scheme is `mcnp`, with `nndc` also supported.

## Remaining limits

The two data workflows remain distinct: multigroup output has `filetype=mgxs`, while continuous-energy neutron output has `filetype=data_neutron`. The bundled multigroup tables cannot recover continuous-energy detail. Continuous-energy parsing and serialization follow OpenMC's capabilities; this repository does not independently reconstruct reaction physics, collapse energy groups, or convert thermal scattering and other particle data. Temperatures merged into a nuclide must come from the same evaluation.

OpenMC neutron HDF5 loading is verified through its Python API. No OpenMC/OpenSn transport runs were executed. Finite-geometry k-effective values and the bundled `pyopensn` driver remain unverified; the checks establish data conversion and loading, not transport benchmark agreement.

MCNP geometry conversion is a surface-based approximation, not a full MCNP input translator. It does not translate cell Boolean expressions, material mappings, transformations, or boundary conditions. Slab transverse extents and infinite-domain proxies need review; sphere and cylinder region builders assume concentric geometry. Mesh checks establish readability and the existence of 3D physical regions, not transport equivalence or a mesh convergence result. STEP import was not tested because no STEP fixture is bundled.

P1 reconstruction from equiprobable angular bins is approximate, higher angular moments are not retained, and invalid angular sections do not receive the same comprehensive validation as the principal cross-section blocks. Input validation is not a complete scientific validation of arbitrary external nuclear data.

Existing example HDF5/mesh files were left untouched. Use newly generated output directories to benefit from the corrected writers. The historical script's imported helper functions remain archival; its CLI and the `multi_group_ace` package use the maintained implementation.
