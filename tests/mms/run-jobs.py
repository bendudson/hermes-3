import os
from xbout import open_boutdataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from perpendicular_laplacian import fstr, astr

def lin_func(x,b,a):
    return b*x + a
    
def collectvar(datasets, var, mesh=0):
    return datasets[mesh][var]

nnbase = 20
ddbase = 0.05
ntest = 3
workdirs = []
# make test for each resolution based on template file
for i in range(0,ntest):
    workdir = f"slab-mms-test-{i}"
    workdirs.append(workdir)
    # create directory 
    if not os.path.isdir(workdir):
        os.system("mkdir "+workdir)
    # copy template
    file = workdir+"/BOUT.inp"
    os.system(f"cp BOUT.inp.template "+file)
    # update with mesh values for test
    nn = nnbase*(i+1)
    dd = ddbase/(i+1)
    with open(file,"a") as file:
        mesh_string = f"""
[mesh]
symmetricGlobalX = false
extrapolate_y = false

nx = {nn}
dx = {dd}
ny = {nn}
dy = {dd}
nz = {nn}
dz = {dd}

g11 = 1.0
g22 = 1.0
g33 = 1.0
g12 = 0.5
g23 = 0.5
g13 = 0.5
x_input = x
y_input = y
z_input = z
a = {astr}
f = {fstr}
"""
        file.write(mesh_string)

    # run job on this input
    print("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")
    os.system("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")

# now analyse the results of the test
# this slice avoids including guard cells in the test
# is this required for periodic geometries?
s = slice(2, -2), slice(None), slice(None)
datasets = []

# open the series of "BOUT.mesh.0.nc" files, 
# saving them in a list `datasets`
for m in range(0,ntest):
    boutmeshpath = workdirs[m]+"/"+f'BOUT.0.nc'
    boutinppath = workdirs[m]+"/"+'BOUT.inp'
    datasets += [ open_boutdataset(boutmeshpath, inputfilepath=boutinppath, keep_yboundaries=False, is_mms_dump=True) ]

#   l2norm = []
#   nylist = []
#   dylist = []
#   for m in range(0,ntest):
    #  numerical = collectvar(datasets, f"result", m)
    #  #print(numerical)
    #  attrs = numerical.attrs
    #  ops, inp = attrs["operator"], attrs["inp"]
    #  #print(ops)
    #  #print(inp)
    #  #print(get_analytical(ops, inp))
    #  analytical = get_analytical(ops, inp)(Rs[m], Zs[m])
    #  #print(analytical)
    #  error_values = (numerical - analytical)[s]
    #  #print(error_values)

    #  thisl2 = np.sqrt(np.mean(error_values**2))
    #  l2norm.append(thisl2)
    #  out[inp][ops].append(thisl2)
    #  nylist.append(Rs[m].shape[1])
    #  # proxy for grid spacing
    #  dylist.append(1.0/Rs[m].shape[1])
