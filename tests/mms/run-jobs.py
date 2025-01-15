import os
from xbout import open_boutdataset
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from perpendicular_laplacian import fstr, astr, div_a_grad_perp_f_str, div_a_grad_perp_f_func

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
g12 = 0.0
g23 = 0.0
g13 = 0.0
x_input = x
y_input = y
z_input = z
a = {astr}
f = {fstr}
expected_result = {div_a_grad_perp_f_str}
"""
        file.write(mesh_string)

    # run job on this input
    print("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")
    os.system("../.././hermes_mms_tests -d "+workdir+" > "+workdir+"/output.txt")

# now analyse the results of the test
# this slice avoids including guard cells in the test
# is this required for periodic geometries?
s = slice(2, -2), slice(None), slice(None)
plot_data = dict()
datasets = []
# open the series of "workdir/BOUT.0.nc" files, 
# saving them in a list `datasets`
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
    analytical = div_a_grad_perp_f_func(xx,yy,zz)
    #print(analytical.values[:,0,0])
    #print(numerical.values[:,0,0])
    #print(expected.values[:,0,0])
    error_values = (numerical - analytical)[s]
    test_error_values = (expected - analytical)[s]
    testl2 = np.sqrt(np.mean(test_error_values**2))
    testl2norm.append(testl2)
    thisl2 = np.sqrt(np.mean(error_values**2))
    l2norm.append(thisl2)
    nylist.append(yy.shape[1])
    # proxy for grid spacing
    dylist.append(1.0/yy.shape[1])
        
# cast lists as numpy arrays for further manipulation
testl2norm = np.array(testl2norm)
print("test analytical error: ",testl2norm)
l2norm = np.array(l2norm)
nylist = np.array(nylist)
dylist = np.array(dylist)

# find linear fit coefficients to test convergence rate
# and construct fit function for plotting
# print(l2norm)
logl2norm = np.log(l2norm)
logdylist = np.log(dylist)
outvals = curve_fit(lin_func,logdylist,logl2norm)
coeffs = outvals[0]
slope = coeffs[0] # also the convergence order
offset = coeffs[1]
logfit = slope*logdylist + offset
fitfunc = np.exp(logfit)

# record results in dictionary and plot
#label = attrs["operator"] + " : f = " + attrs["inp"]
label = "FV::Div_a_Grad_perp(a, f)"
plot_data[label] = [dylist, l2norm, fitfunc, slope, offset]
for key, variable_set in plot_data.items():
    (xaxis, yaxis, fit, slope, offset) = variable_set
    plt.figure()
    plt.plot(xaxis, yaxis, "x-", label="$\\epsilon(\\mathcal{L}\\ast f)$: "+label)
    plt.plot(xaxis, yaxis[0]*(xaxis/xaxis[0])**2, "x-", label="$\\propto \\Delta^2$")
    plt.plot(xaxis, fit, "x-", label="$e^{{{:.2f}}}\\Delta^{{{:.2f}}}$".format(offset,slope))
    plt.xlabel("$\\Delta = 1/N_y$")
    plt.title(key)
    plt.legend()
    plt.gca().set_yscale("log")
    plt.gca().set_xscale("log")
plt.show()
