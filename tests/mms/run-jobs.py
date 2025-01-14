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
# make test for each resolution based on template file
for i in range(0,ntest):
    workdir = f"slab-mms-test-{i}"
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
