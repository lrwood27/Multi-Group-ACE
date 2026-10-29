# Re-export the small public API from the package root so user code can
# import from "multi_group_ace" directly instead of each internal module.
from .core import MGData, KB_EV_PER_K, MEV_TO_EV
from .json_loader import load_from_json
from .converter import build_ace
from .ace_io import read_mg_ace, write_ace
from .hdf5_io import write_mgxs_hdf5_isotropic
from .gmsh_geometry import mcnp_to_gmsh_mesh, make_mesh_for_problem
from .continuous_energy import convert_ace_to_openmc_hdf5

__version__ = "0.2.0"

# Limit wildcard imports to the intended public surface.
__all__ = [
    "MGData",
    "KB_EV_PER_K",
    "MEV_TO_EV",
    "load_from_json",
    "build_ace",
    "read_mg_ace",
    "write_ace",
    "write_mgxs_hdf5_isotropic",
    "mcnp_to_gmsh_mesh",
    "make_mesh_for_problem",
    "convert_ace_to_openmc_hdf5",
]
