import os
from pyopensn.context import UseColor
from pyopensn.mesh import FromFileMeshGenerator, PETScGraphPartitioner
from pyopensn.xs import MultiGroupXS
from pyopensn.aquad import GLCProductQuadrature3DXYZ
from pyopensn.solver import DiscreteOrdinatesProblem, PowerIterationKEigenSolver
from pyopensn.fieldfunc import FieldFunctionGridBased

UseColor(False)

HERE = os.path.dirname(os.path.abspath(__file__))
MESH_FILE = os.path.join(HERE, "geom.msh")
XS_FILE   = os.path.join(HERE, "mgxs.h5")
TEMP_K    = 294.0

DATASET_NAME = "99902.02m"

# --- Mesh: read the gmsh file ---
meshgen = FromFileMeshGenerator(
    filename=MESH_FILE,
    partitioner=PETScGraphPartitioner(type="parmetis")  # good for MPI runs
)
grid = meshgen.Execute()
grid.ExportToPVTU("mesh")  # view mesh in ParaView

# If your domain is just the sphere, all boundary faces are the sphere surface.
# Give them ONE boundary name so we can apply vacuum everywhere on the outer surface.
grid.SetUniformBoundaryID("vac")  # name all boundary faces "vac"

# See what material block IDs exist (these come from Gmsh "physical" region IDs)
block_vols = grid.ComputeVolumePerBlockID()
print("Block volumes:", block_vols)

# If you truly have only the Pu region, you'll see one key here.
# Otherwise, you'll need to map each block_id to the correct XS.
PU_BLOCK_ID = sorted(block_vols.keys())[0]

# --- Cross sections: load OpenMC MGXS ---
xs_pu = MultiGroupXS()
xs_pu.LoadFromOpenMC(XS_FILE, DATASET_NAME, TEMP_K)

try:
    ng = int(xs_pu.num_groups)
except Exception:
    ng = 1  # fallback for your 1-group library
print("ng =", ng)

import re

def make_glc_xyz(polar=8, azim=16, scat=0):
    # Defaults keyed by the *possible* names OpenSn might require.
    # If OpenSn asks for a specific key, we fill it from this map.
    defaults = {
        "scattering_order": scat,

        # polar resolution
        "polar_angles": polar,
        "polar_order": polar,
        "num_polar": polar,
        "n_polar": polar,

        # azimuthal resolution
        "azimuthal_angles": azim,
        "azimuthal_order": azim,
        "num_azimuthal": azim,
        "n_azimuthal": azim,

        # sometimes bindings call one of these quadrature_order
        "quadrature_order": polar,
    }

    kwargs = {}
    while True:
        try:
            q = GLCProductQuadrature3DXYZ(**kwargs)
            print("GLCProductQuadrature3DXYZ kwargs used:", kwargs)
            return q
        except RuntimeError as e:
            m = re.search(r'Key "([^"]+)" must be provided', str(e))
            if not m:
                raise
            key = m.group(1)
            if key not in defaults:
                print("GLCProductQuadrature3DXYZ init doc:", getattr(GLCProductQuadrature3DXYZ.__init__, "__doc__", None))
                raise RuntimeError(f"OpenSn requires key '{key}', but I don't have a default for it yet.")
            kwargs[key] = defaults[key]

quad = make_glc_xyz(polar=8, azim=16, scat=xs_pu.scattering_order)

try:
    ng = int(xs_pu.num_groups)
except Exception:
    ng = 1
print("ng =", ng)

quad = make_glc_xyz(polar=8, azim=16, scat=int(xs_pu.scattering_order))

groupsets = [{
    "groups_from_to": [0, ng - 1],
    "angular_quadrature": quad,
    "inner_linear_method": "petsc_gmres",
    "l_abs_tol": 1.0e-6,
    "l_max_its": 200,
    "gmres_restart_interval": 30,
    "angle_aggregation_type": "single",
    "angle_aggregation_num_subsets": 1,
}]

# Map mesh block IDs to XS objects
xs_map = [{"block_ids": [PU_BLOCK_ID], "xs": xs_pu}]

# Some OpenSn builds require xs_map at construction time
problem = DiscreteOrdinatesProblem(mesh=grid, num_groups=ng, groupsets=groupsets, xs_map=xs_map)

# Vacuum boundary condition on boundary name vac
problem.SetBoundaryOptions(boundary_conditions=[{"name": "vac", "type": "vacuum"}])


# --- k-eigen solver (power iteration) ---
k_solver = PowerIterationKEigenSolver(problem=problem, max_iters=200, k_tol=1.0e-8)
k_solver.Initialize()
k_solver.Execute()

print("k-effective:", k_solver.GetEigenvalue())

# --- Export scalar fluxes (one field function per energy group) ---
phi_g = problem.GetScalarFieldFunctionList(only_scalar_flux=True)
FieldFunctionGridBased.ExportMultipleToPVTU(phi_g, "phi")

print("Wrote mesh.pvtu and phi_*.pvtu for ParaView.")
