# Copilot Instructions for Multi-Group-ACE

## Project Overview
**Multi-Group-ACE** is a modular Python package for multigroup workflows and conversion of existing continuous-energy neutron ACE libraries into OpenMC HDF5. The maintained entry point is `python -m multi_group_ace.cli`; executing the historical `simple_ace_mg_full.py` delegates to the same CLI. Its old helper definitions remain archival.

The code snippets below are multigroup teaching sketches from the original refactor, not complete current implementations to paste into the package. The source modules, README, and technical reference define the current interfaces. Continuous-energy conversion uses OpenMC objects directly and does not pass through `MGData`.

## Research Meeting 1 Requirements

### What to Remove
- ❌ **`@dataclass` decorators**: Replace with plain `class` + `__init__`
- ❌ **`from typing import ...`**: Remove ALL type hints (`List[float]`, `Optional[str]`, `Tuple[...]`)
- ❌ **Inline type comments**: Convert to block comments above code
- ❌ **Multiple JSON key fallbacks**: One canonical name per field (no `sigt OR sigma_t OR Sigma_t`)
- ❌ **Excessive function parameters**: Keep essential args; split variants into separate functions

### What to Add
- ✅ **Plain classes** with `__init__` methods
- ✅ **Block comments** (standard Python style)
- ✅ **`if __name__ == "__main__"`** entry point
- ✅ **`__init__.py`** with module exports

## Current Architecture (9 Modules)

```
multi_group_ace/
├── __init__.py                 # Package exports
├── cli.py                      # Command-line interface + main()
├── core.py                     # MGData class, constants
├── json_loader.py              # Load config from JSON (strict keys)
├── converter.py                # build_ace() - ACE table generation
├── ace_io.py                   # read_mg_ace(), write_ace()
├── hdf5_io.py                  # write_mgxs_hdf5_isotropic()
├── continuous_energy.py        # OpenMC neutron ACE -> HDF5 + XML
└── gmsh_geometry.py            # MCNP geometry → Gmsh mesh
```

## File Breakdown (What Each Does)

### 1. `core.py` — Data Container & Constants
**Purpose**: Define MGData class and physical constants  
**Size**: ~40 lines  
**Key Concepts**:
- MGData: Plain class holding cross-section data
- Constants: KB_EV_PER_K, MEV_TO_EV (temperature conversion)
- Energy grids always descending (MeV)

**What's Inside**:
```python
# Physical constants
KB_EV_PER_K = 8.617333262e-5
MEV_TO_EV = 1.0e6

# Plain class (no @dataclass)
class MGData:
    def __init__(self, zaid, NG, erg, sigt, sigc, sigf, nu, chi, sigs, sig1=None):
        self.zaid = zaid        # Material ID (e.g., "99902.01m")
        self.NG = NG            # Number of energy groups
        self.erg = erg          # NG+1 energy boundaries (descending, MeV)
        self.sigt = sigt        # Total cross-section (NG values)
        self.sigc = sigc        # Capture cross-section
        self.sigf = sigf        # Fission cross-section
        self.nu = nu            # Neutrons per fission (NG values)
        self.chi = chi          # Fission spectrum (NG values)
        self.sigs = sigs        # Scattering matrix (NG*NG, row-major)
        self.sig1 = sig1 or []  # P1 anisotropic (optional, often empty)
```

**Explain to Professor**: "MGData is the shared container for multigroup JSON, ACE, and MGXS workflows. Continuous-energy ACE conversion instead uses OpenMC IncidentNeutron objects to retain its energy grids and reaction data."

---

### 2. `json_loader.py` — Configuration Input
**Purpose**: Load nuclear data config from JSON file  
**Size**: ~60 lines  
**Key Concepts**:
- ONE canonical key name per field (strict validation)
- Early error detection if keys missing
- Micro/macro cross-section toggle

**What's Inside**:
```python
import json
from core import MGData

def load_from_json(path):
    # Load multigroup config from JSON file
    # Expects canonical keys: zaid, NG, erg, sigt, sigc, sigf, nu, chi, sigs
    # Raises KeyError immediately if required key missing (no fallbacks)
    
    with open(path) as f:
        data = json.load(f)
    
    # Validate required keys
    required_keys = ["zaid", "NG", "erg", "sigt", "sigc", "sigf", "nu", "chi", "sigs"]
    for key in required_keys:
        if key not in data:
            raise KeyError(f"Missing required key in JSON: {key}")
    
    # Optional: check if micro or macro cross-sections
    is_micro = data.get("is_micro", True)
    
    # Create MGData object
    mg = MGData(
        zaid=data["zaid"],
        NG=data["NG"],
        erg=data["erg"],
        sigt=data["sigt"],
        sigc=data["sigc"],
        sigf=data["sigf"],
        nu=data["nu"],
        chi=data["chi"],
        sigs=data["sigs"],
        sig1=data.get("sig1", None)
    )
    
    return mg, is_micro
```

**Explain to Professor**: "This function reads a JSON file and converts it into our MGData object. It's strict—if a required field is missing, it fails immediately instead of trying multiple spellings."

---

### 3. `converter.py` — Physics Core (ACE Generation)
**Purpose**: Build ACE table arrays from MGData  
**Size**: ~150 lines  
**Key Concepts**:
- ACE format: three arrays (nxs, jxs, xss)
- 1-based indexing (ACE quirk; confusing but standard)
- Equiprobable µ-bin calculation
- Energy grid setup

**What's Inside**:
```python
from core import MGData

def get_equiprob_bins(nbins, N, A):
    # Generate equiprobable µ (scattering angle cosine) bins
    # Input: nbins=number of bins, N=number of scattering reactions, A=cumulative PDF
    # Output: bin boundaries for µ
    # (Port of Perl get_equiprob_bins from original script)
    
    bins = []
    for i in range(nbins + 1):
        mu = -1.0 + 2.0 * i / nbins
        bins.append(mu)
    return bins

def build_ace(mg):
    # Build ACE cross-section tables from MGData
    # Returns: (nxs, jxs, xss) lists
    # nxs: header integers (material info)
    # jxs: pointers to data sections in xss
    # xss: actual cross-section values
    
    NG = mg.NG
    nxs = [0] * 16
    jxs = [0] * 32
    xss = []
    
    # NXS array setup
    nxs[0] = NG                      # Number of energy groups
    nxs[1] = 1                       # Always 1 for ACE
    
    # Build energy grid (descending)
    erg_section = list(reversed(mg.erg))  # Convert to ascending for ACE
    xss.extend(erg_section)
    
    # Build cross-section sections (sigt, sigc, sigf, etc.)
    jxs[1] = len(xss) + 1            # Pointer to total XS (1-based)
    xss.extend(mg.sigt)
    
    jxs[2] = len(xss) + 1            # Pointer to capture XS
    xss.extend(mg.sigc)
    
    jxs[3] = len(xss) + 1            # Pointer to fission XS
    xss.extend(mg.sigf)
    
    # Add nu, chi, scattering, etc.
    jxs[4] = len(xss) + 1
    xss.extend(mg.nu)
    
    jxs[5] = len(xss) + 1
    xss.extend(mg.chi)
    
    jxs[6] = len(xss) + 1
    xss.extend(mg.sigs)
    
    return nxs, jxs, xss
```

**Explain to Professor**: "This is the core physics module. It takes the MGData structure and converts it into ACE format arrays. ACE uses 1-based indexing (pointer arrays) which is why the code looks confusing—it's just tracking where each data section starts in the xss array."

---

### 4. `ace_io.py` — ACE File I/O
**Purpose**: Read and write ACE format files  
**Size**: ~120 lines  
**Key Concepts**:
- ACE binary format (32-bit floats, 1-based indexing)
- Text format support (optional)
- Round-trip: file → MGData → file

**What's Inside**:
```python
import struct
from core import MGData

def read_mg_ace(ace_file):
    # Read ACE cross-section file and return MGData object
    # Supports simple multigroup ACE files (ASCII or binary)
    
    mg = None
    with open(ace_file, 'r') as f:
        lines = f.readlines()
    
    # Parse ACE header and tables (simplified for multigroup)
    # Extract NG, erg, sigt, sigc, sigf, nu, chi, sigs
    # Build MGData from parsed data
    
    return mg

def write_ace(nxs, jxs, xss, output_file):
    # Write ACE arrays to file (ASCII format for readability)
    
    with open(output_file, 'w') as f:
        # Write NXS array
        f.write("NXS\n")
        for val in nxs:
            f.write(f"{val}\n")
        
        # Write JXS array
        f.write("JXS\n")
        for val in jxs:
            f.write(f"{val}\n")
        
        # Write XSS array
        f.write("XSS\n")
        for val in xss:
            f.write(f"{val}\n")
```

**Explain to Professor**: "These functions handle file I/O for ACE format. Read converts a file on disk back into MGData; write saves the ACE tables to disk. This lets us chain workflows (JSON → ACE file → HDF5)."

---

### 5. `hdf5_io.py` — Multigroup MGXS Output
**Purpose**: Write mgxs.h5 files compatible with OpenMC  
**Size**: ~100 lines  
**Key Concepts**:
- HDF5 hierarchical format
- Temperature in eV (convert from MeV)
- Isotropic scattering (no P1 for now)
- OpenMC mgxs schema

**What's Inside**:
```python
try:
    import h5py
except ImportError:
    h5py = None

from core import MGData, KB_EV_PER_K, MEV_TO_EV

def write_mgxs_hdf5_isotropic(mg, h5_path, library_name=None):
    # Write OpenMC-compatible mgxs.h5 file from MGData
    # Isotropic scattering (no P1)
    
    if h5py is None:
        raise ImportError("h5py required for HDF5 output: pip install h5py")
    
    with h5py.File(h5_path, 'w') as f:
        # Create library metadata
        lib = f.create_group('library')
        lib.attrs['name'] = library_name or 'multigroup'
        lib.attrs['version'] = '1.0'
        
        # Create material group
        mat = f.create_group(f'material_{mg.zaid}')
        mat.attrs['zaid'] = mg.zaid
        mat.attrs['n_groups'] = mg.NG
        
        # Store cross-sections
        mat.create_dataset('energy_groups', data=mg.erg)
        mat.create_dataset('total_xs', data=mg.sigt)
        mat.create_dataset('capture_xs', data=mg.sigc)
        mat.create_dataset('fission_xs', data=mg.sigf)
        mat.create_dataset('nu', data=mg.nu)
        mat.create_dataset('chi', data=mg.chi)
        mat.create_dataset('scatter_xs', data=mg.sigs)
        
        # Temperature (convert MeV → eV, compute K)
        if hasattr(mg, 'tmp'):
            tmp_ev = mg.tmp * MEV_TO_EV
            tmp_k = tmp_ev / KB_EV_PER_K
            mat.attrs['temperature_K'] = tmp_k
```

**Explain to Professor**: "This module writes OpenMC/OpenSn multigroup MGXS from MGData. It reorders scattering arrays and converts ACE kT from MeV to eV and Kelvin labels, retaining each material's temperature. P1 uses order 1 while the incoming-direction representation remains isotropic. Continuous-energy neutron HDF5 uses a separate writer."

---

### 6. `continuous_energy.py` — OpenMC Continuous-Energy Neutron Output

**Purpose**: Convert existing continuous-energy neutron ACE evaluations into native OpenMC HDF5 plus an XML catalog. OpenMC 0.16.0 is optional and loaded only for this workflow.

```text
Neutron ACE file(s), every table enumerated
    -> OpenMC IncidentNeutron objects
    -> compatible temperatures grouped by nuclide
    -> <GNDS nuclide>.h5 + cross_sections.xml
```

The converter detects table types from contents and delegates supported ASCII/binary neutron parsing and HDF5 serialization to OpenMC. It does not reconstruct continuous-energy information from multigroup values. The bundled `*.m` examples remain multigroup data; thermal scattering and other particle tables are unsupported.

```python
from multi_group_ace import convert_ace_to_openmc_hdf5

hdf5_paths, xml_path = convert_ace_to_openmc_hdf5(
    ["He3.ace", "U235.ace"], "openmc_data",
    metastable_scheme="mcnp", overwrite=False,
)
```

The result contains absolute paths. Each nuclide gets one `filetype=data_neutron` file with compatible temperatures combined without averaging cross sections. Merged temperatures must originate from the same evaluation. Duplicate rounded temperature labels, incompatible reaction sets, and incompatible atomic weight ratios are rejected.

Use `--ace FILE [FILE ...] --outdir openmc_data` for multiple neutron files; the default output directory is `h5_out`. A custom `--h5` requires exactly one resulting nuclide and a `.h5` suffix, cannot be combined with an explicit `--outdir`, and places XML in the output file's parent. `--metastable-scheme` accepts `mcnp` (default) or `nndc`. Existing HDF5/XML destinations require `--overwrite`; the replacement catalog contains only this invocation's outputs. `--library-name` remains multigroup-only.

**Explain to Professor**: "Continuous-energy ACE carries reaction and distribution information absent from MGData. We use OpenMC's data model and exporter, enumerate all source tables, validate temperature merges, and generate the catalog OpenMC uses to locate the neutron files."

---

### 7. `gmsh_geometry.py` — MCNP → Mesh Conversion
**Purpose**: Parse MCNP geometry cards and generate Gmsh meshes  
**Size**: ~180 lines  
**Key Concepts**:
- MCNP surface cards (PX, PY, CY, SP, etc.)
- Gmsh Python API or file generation
- Region definition from surfaces

**What's Inside**:
```python
import re
import subprocess

class Plane:
    # Represents a plane surface (PX, PY, PZ, etc.)
    def __init__(self, name, normal, position):
        self.name = name
        self.normal = normal        # (0,1,0) for YZ plane, etc.
        self.position = position    # Offset distance
    
    def to_gmsh_geo(self):
        # Return Gmsh .geo file syntax for this plane
        pass

class Cylinder:
    # Represents a cylindrical surface (CX, CY, CZ, etc.)
    def __init__(self, name, axis, radius):
        self.name = name
        self.axis = axis            # 'x', 'y', 'z'
        self.radius = radius
    
    def to_gmsh_geo(self):
        # Return Gmsh syntax
        pass

def parse_mcnp_surfaces(geo_file):
    # Parse MCNP geometry file (.inp or .geo)
    # Extract surface definitions (PX, PY, CY, SP, etc.)
    # Return list of Plane/Cylinder/Sphere objects
    
    surfaces = []
    with open(geo_file, 'r') as f:
        for line in f:
            # Parse "PX 1.5" → Plane('PX', normal, position)
            match = re.match(r'(P[XYZ]|C[XYZ]|SP)\s+([\d.eE+-]+)', line)
            if match:
                surf_type = match.group(1)
                value = float(match.group(2))
                surfaces.append((surf_type, value))
    
    return surfaces

def mcnp_to_gmsh_mesh(geo_file, output_msh):
    # Parse MCNP geometry, generate Gmsh .geo file, run gmsh
    # Output: geom.msh file
    
    surfaces = parse_mcnp_surfaces(geo_file)
    
    # Write temporary .geo file
    with open('_temp.geo', 'w') as f:
        f.write("SetFactory('OpenCASCADE');\n")
        for surf_type, value in surfaces:
            if surf_type == 'PX':
                f.write(f"Box({value}, -10, -10, 0.1, 20, 20);\n")
            # ... handle other surface types
    
    # Run Gmsh
    subprocess.run(['gmsh', '_temp.geo', '-3', '-o', output_msh])
```

**Explain to Professor**: "This module bridges MCNP geometry format and Gmsh mesh generation. It parses MCNP surface cards (planes, cylinders) and converts them into Gmsh syntax. Then it calls the Gmsh command-line tool to generate a 3D mesh file."

---

### 8. `cli.py` — Command-Line Interface
**Purpose**: Argument parsing and workflow orchestration  
**Size**: ~100 lines  
**Key Concepts**:
- `argparse` for CLI arguments
- Three top-level modes route four data workflows: JSON to multigroup ACE, multigroup ACE to MGXS, neutron ACE to continuous-energy HDF5, and multigroup verification scanning
- `if __name__ == "__main__"` entry point

**What's Inside**:
```python
import argparse
import sys
from json_loader import load_from_json
from converter import build_ace
from ace_io import write_ace, read_mg_ace
from hdf5_io import write_mgxs_hdf5_isotropic
from gmsh_geometry import mcnp_to_gmsh_mesh

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Multi-Group ACE processor'
    )
    
    # Workflow 1: JSON to ACE
    parser.add_argument('--config', help='Load from JSON config file')
    
    # Workflow 2: ACE to HDF5
    parser.add_argument('--ace', help='Input ACE file')
    parser.add_argument('--h5', help='Output HDF5 file')
    
    # Workflow 3: Verification tree scan
    parser.add_argument('--verification-root', help='Root directory of problems')
    parser.add_argument('--outdir', help='Output directory for HDF5 files')
    parser.add_argument('--make-mesh', action='store_true', help='Generate Gmsh meshes')
    parser.add_argument('--gmsh', default='gmsh', help='Gmsh executable path')
    
    args = parser.parse_args()
    
    # Execute appropriate workflow
    if args.config:
        mg, is_micro = load_from_json(args.config)
        nxs, jxs, xss = build_ace(mg)
        output_file = args.config.replace('.json', '.ace')
        write_ace(nxs, jxs, xss, output_file)
        print(f"Wrote ACE to {output_file}")
    
    elif args.ace and args.h5:
        mg = read_mg_ace(args.ace)
        write_mgxs_hdf5_isotropic(mg, args.h5)
        print(f"Wrote HDF5 to {args.h5}")
    
    elif args.verification_root:
        # Tree scan workflow
        pass

if __name__ == "__main__":
    main()
```

**Explain to Professor**: "The CLI validates three mutually exclusive modes. ACE mode detects contents and chooses the multigroup or continuous-energy writer. Recursive scanning remains multigroup-only; neutron ACE mode processes multiple files and temperatures. The simplified code sketch above shows only the multigroup branches."

---

### 9. `__init__.py` — Package Interface
**Purpose**: Export public API  
**Size**: ~10 lines  
**What's Inside**:
```python
# multi_group_ace/__init__.py

from core import MGData, KB_EV_PER_K, MEV_TO_EV
from json_loader import load_from_json
from converter import build_ace
from ace_io import read_mg_ace, write_ace
from hdf5_io import write_mgxs_hdf5_isotropic
from gmsh_geometry import mcnp_to_gmsh_mesh

__version__ = '0.1.0'
__all__ = [
    'MGData',
    'KB_EV_PER_K',
    'MEV_TO_EV',
    'load_from_json',
    'build_ace',
    'read_mg_ace',
    'write_ace',
    'write_mgxs_hdf5_isotropic',
    'mcnp_to_gmsh_mesh',
]
```

**Explain to Professor**: "This file defines the public API. Users can `from multi_group_ace import MGData` instead of `from multi_group_ace.core import MGData`. It's cleaner and abstracts the internal structure."

---

## Current Continuous-Energy Verification

Use the optional API in an isolated Windows environment:

```powershell
python -m venv --system-site-packages .venv-ce
.venv-ce\Scripts\python -m pip install -r requirements-ce.txt
.venv-ce\Scripts\python tests/fetch_ce_data.py
$env:MULTI_GROUP_ACE_REQUIRE_CE_TESTS = "1"
.venv-ce\Scripts\python -m unittest discover -s tests -v
```

The pinned OpenMC 0.16.0 API requires Python 3.12 or newer and Git/network access for its source installation; verification used Python 3.13.15. No transport executable is needed. The downloader verifies authentic He3 and U235 ENDF/B-VII.1 fixtures from pinned NJOY ACEtk resources. Downloads, tools, and isolated environments remain ignored. Real integration tests compare outputs with direct OpenMC ACE loading and reopen HDF5/XML catalogs. They also cover binary ACE, concatenated tables, multiple files, and temperature merging. Modified-temperature inputs are explicitly synthetic plumbing tests.

Retain all 24 original multigroup checks. Dependency-light runs may skip optional checks; the complete run before committing must have no skips. Record actual results in `doc/repository_audit.md`. Data-conversion checks do not establish finite-geometry transport benchmark agreement.

The complete verification passed 51 tests with zero skips or failures: 24 original workflow checks, 14 continuous-energy unit checks, and 13 real OpenMC integration checks. The required-CE setting makes missing OpenMC or fixtures an error.

## How to Explain the Current Package

### 30-Second Version
*"The maintained package has nine focused modules. It handles multigroup data through MGData and uses OpenMC's Python API for continuous-energy neutron ACE conversion. The two physical representations stay separate while sharing one CLI and repeatable verification."*

### 2-Minute Version
*"The multigroup path loads JSON into MGData, builds or reads ACE arrays, and writes MGXS HDF5. The continuous-energy path reads every neutron ACE table through OpenMC, checks compatible temperatures, and writes one neutron HDF5 file per nuclide plus cross_sections.xml. The mesh module builds limited surface-based geometry, and the CLI routes each request. These conversions do not execute a transport solver or reconstruct information missing from grouped data."*

### Visual Summary
```
JSON config
    ↓ (json_loader.py)
MGData object
    ↓ (converter.py) ← Physics core: build_ace()
[nxs, jxs, xss] arrays
    ↓ (ace_io.py)
ACE file on disk
    ↓ (read_mg_ace)
MGData
    ↓ (hdf5_io.py) ← Temperature conversion (MeV → eV)
mgxs.h5 (OpenMC format)

Optional: geometry_file → gmsh_geometry.py → geom.msh
```

## Historical Refactoring Checklist

These extraction tasks describe the original refactor, which is already present in the maintained package. They remain as teaching context. Use the current verification procedure below after changes.

### Phase 1: Create Modules (1–2 hours)
- [ ] Create `multi_group_ace/` directory
- [ ] Extract `core.py` (MGData class + constants)
- [ ] Extract `json_loader.py` (load_from_json)
- [ ] Extract `converter.py` (build_ace, get_equiprob_bins)
- [ ] Extract `ace_io.py` (read_mg_ace, write_ace)
- [ ] Extract `hdf5_io.py` (write_mgxs_hdf5_isotropic)
- [ ] Extract `gmsh_geometry.py` (MCNP parsing, mesh gen)
- [ ] Extract `cli.py` (argparse, main)
- [ ] Create `__init__.py`

### Phase 2: Remove Type Hints & Dataclasses (30 min)
- [ ] Replace `@dataclass` with `class MGData.__init__`
- [ ] Remove all `from typing import ...` lines
- [ ] Remove all `: Type` annotations from functions
- [ ] Replace inline `# type:` comments with block comments

### Phase 3: Test (30 min)
- [ ] `python -m multi_group_ace.cli --config json examples/ex3.json`
- [ ] `python -m multi_group_ace.cli --ace test.ace --h5 test.h5`
- [ ] Verify outputs match original script

---

*Historical refactoring checklist recorded January 28, 2026. Current usage and verification updated October 29, 2026.*
