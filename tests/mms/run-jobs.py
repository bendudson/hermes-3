import os
from xbout import open_boutdataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from perpendicular_laplacian import fstr, astr, div_a_grad_perp_f_str, div_a_grad_perp_f_func
from perpendicular_laplacian import g11_str, g12_str, g13_str, g22_str, g23_str, g33_str

def lin_func(x,b,a):
    return b*x + a
    
def collectvar(datasets, var, mesh=0):
    return datasets[mesh][var]

# This test script below works on BOUT.0.nc output with a single 
# operator tested. This script can be generalised to work on multiple operators.

nnbase = 20 # base number of grid points
ddbase = 0.05 # base grid spacing
ntest = 3 # number of grids tested
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
    nn = nnbase*(i+1) # number of points in each grid
    dd = 2.0*np.pi/nn # y z grid spacing, z, y on [0, 2pi]
    ddx = 1.0/(nn-4) # x grid spacing to account for guard cells, x on [0,1] 
    # symmetricGlobalX = true ensures that x = 0 and x = 1 sits
    # halfway between the last true grid point and the first guard point.
    with open(file,"a") as file:
        mesh_string = f"""
[mesh]
symmetricGlobalX = true
extrapolate_y = false
extrapolate_x= false
extrapolate_z= false

nx = {nn}
dx = {ddx}
ny = {nn}
dy = {dd}
nz = {nn}
dz = {dd}

g11 = {g11_str}
g22 = {g22_str}
g33 = {g33_str}
g12 = {g12_str}
g23 = {g23_str}
g13 = {g13_str}
x_input = x
y_input = y
z_input = z
a = {astr}
f = {fstr}
expected_result = {div_a_grad_perp_f_str}
"""
        file.write(mesh_string.replace("**","^"))

    # run job on this input
    print("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")
    os.system("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")

# now analyse the results of the test
# this slice avoids including guard cells in the test
# need guard cells in x (assume 2 here) and guard cells in y (assume 1)
# no guard cells in z
s = slice(2, -2), slice(1, -1), slice(None)

# a dictionary of plot data, filled later on
plot_data = dict()

# open the series of "workdir/BOUT.0.nc" files, 
# saving them in a list `datasets`
datasets = []
for workdir in workdirs:
    boutmeshpath = workdir+"/"+f'BOUT.0.nc'
    boutinppath = workdir+"/"+'BOUT.inp'
    datasets.append(open_boutdataset(boutmeshpath, inputfilepath=boutinppath, keep_yboundaries=False))

testl2norm = []
l2norm = []
nylist = []
dylist = []
for m in range(0,ntest):
    numerical = collectvar(datasets, "result", m)
    expected = collectvar(datasets, "expected_result", m)
    xx = collectvar(datasets, "x_input", m)
    yy = collectvar(datasets, "y_input", m)
    zz = collectvar(datasets, "z_input", m)
    ff = collectvar(datasets, "f", m)
    aa = collectvar(datasets, "a", m)
    analytical = div_a_grad_perp_f_func(xx,yy,zz)
    #print(analytical.values[:,5,5])
    #print("num:",numerical.values[:,5,5])
    #print("exp:",expected.values[:,5,5])
    #print("xx:",xx.values[:,5,5])
    #print("a:",aa.values[:,5,5])
    #print("f:",ff.values[:,5,5])
    #print("num:",numerical.values[5,5,:])
    #print("exp:",expected.values[5,5,:])
    #print("zz:",zz.values[5,5,:])
    #print("a:",aa.values[5,5,:])
    #print("f:",ff.values[5,5,:])
    #print("num:",numerical.values[5,:,5])
    #print("exp:",expected.values[5,:,5])
    #print("yy:",yy.values[5,:,5])
    #print("a:",aa.values[5,:,5])
    #print("f:",ff.values[5,:,5])
    error_values = (numerical - analytical)[s]
    #print("err:",error_values.values[:,5,5])
    test_error_values = (expected - analytical)[s]
    testl2 = np.sqrt(np.mean(test_error_values**2))
    testl2norm.append(testl2)
    thisl2 = np.sqrt(np.mean(error_values**2))
    #thisl2 = np.sqrt(np.mean(error_values[:,5,5]**2))
    print(thisl2)
    l2norm.append(thisl2)
    nylist.append(yy.shape[1])
    # proxy for grid spacing
    dylist.append(1.0/yy.shape[1])
        
# cast lists as numpy arrays for further manipulation
testl2norm = np.array(testl2norm)
print("test analytical error: ",testl2norm)
l2norm = np.array(l2norm)
print("test error: ",l2norm)
nylist = np.array(nylist)
dylist = np.array(dylist)

# find linear fit coefficients to test convergence rate
# and construct fit function for plotting
try:
    logl2norm = np.log(l2norm)
    logdylist = np.log(dylist)
    outvals = curve_fit(lin_func,logdylist,logl2norm)
    coeffs = outvals[0]
    slope = coeffs[0] # the convergence order
    offset = coeffs[1]
    logfit = slope*logdylist + offset
    fitfunc = np.exp(logfit)
except ValueError:
    print("Infs/Nans encountered in fit, skipping")
    slope = None
    offset = None
    fitfunc = None

# record results in dictionary and plot
#label = attrs["operator"] + " : f = " + attrs["inp"]
label = "FV::Div_a_Grad_perp(a, f)"
plot_data[label] = [dylist, l2norm, fitfunc, slope, offset]

for key, variable_set in plot_data.items():
    (xaxis, yaxis, fit, slope, offset) = variable_set
    plt.figure()
    plt.plot(xaxis, yaxis, "x-", label="$\\epsilon(\\mathcal{L}\\ast f)$: "+label)
    plt.plot(xaxis, yaxis[0]*(xaxis/xaxis[0])**2, "x-", label="$\\propto \\Delta^2$")
    if not fit is None:
        plt.plot(xaxis, fit, "x-", label="$e^{{{:.2f}}}\\Delta^{{{:.2f}}}$".format(offset,slope))
    plt.xlabel("$\\Delta = 1/N_y$")
    plt.title(key)
    plt.legend()
    if not fit is None:
        plt.gca().set_yscale("log")
        plt.gca().set_xscale("log")
    else:
        print("l2 error: ",yaxis)
plt.show()

