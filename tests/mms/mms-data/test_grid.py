from hypnotoad.cases.circular import CircularEquilibrium
from hypnotoad.core.mesh import BoutMesh
import numpy as np
import copy

R0 = 2.3
B0 = 3.2
r_inner = 0.1
r_outer = 1.4
q = 4.1
settings_default = {
    "R0": R0,
    "B0": B0,
    "poloidal_spacing_method": "linear",
    "q_coefficients": [q],
    "r_inner": r_inner,
    "r_outer": r_outer,
    "refine_methods": "line",
    "refine_width": 1.0e-3,
    "single_region" : True,
    "nx" : 20,
    "ny" : 20,
    "orthogonal" : True,
}

nn_list = [20,40,80] 

for j, nn in enumerate(nn_list):
    # modify the settings for these resolutions
    settings = copy.deepcopy(settings_default)
    settings["nx"] = nn
    settings["ny"] = nn
    # generate the grid file
    equilib = CircularEquilibrium(settings, nonorthogonal_settings=settings)
    mesh = BoutMesh(equilib, settings)
    mesh.geometry()
    mesh.writeGridfile(f"circ_grid_{j}.nc")
