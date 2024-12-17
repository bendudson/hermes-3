from xbout import open_boutdataset
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import collections
import glob
import pickle


def get_convergence(method):
    return 2


expected_failures = {"Div_a_Grad_perp_nonorthog(1, f)", "FV::Div_a_Grad_perp(1, f)"}

success = True
failed = collections.defaultdict(list)
failed2 = set()


def collectvar(datasets, var, mesh=0):
    return datasets[mesh][var]

def r_func(R, Z):
    R0 = 5
    theta = np.arctan2(Z, R - R0)
    beta = R0 * np.cos(2 * theta)
    gamma = 2 * (np.cos(theta) ** 3) - np.cos(theta)
    r = np.sqrt((R - R0) ** 2 + Z**2)
    return (beta + 2 * r * gamma) / (r * beta + r**2 * gamma)


def sin_theta(R, Z):
    R0 = 5
    theta = np.arctan2(Z, R - R0)
    r = np.sqrt((R - R0) ** 2 + Z**2)
    st = np.sin(theta)
    ct = np.cos(theta)
    c2t = np.cos(theta * 2)
    s2t = np.sin(theta * 2)
    A = st - 6 * c2t * st - 2 * R0 / r * s2t
    B = r * R0 * c2t + r**2 * (2 * ct**3 - ct)
    return A / B * ct - st / r**2


ana_default = {
    "R": lambda R, Z: 1 / R,
    "R²": lambda R, Z: R * 0 + 4,
    "sin(R)": lambda R, Z: np.cos(R) / R - np.sin(R),
    "sin(10*R)": lambda R, Z: 10 * np.cos(10 * R) / R - 1e2 * np.sin(10 * R),
    "sin(100*R)": lambda R, Z: 100 * np.cos(100 * R) / R - 1e4 * np.sin(100 * R),
    "sin(1000*R)": lambda R, Z: 1000 * np.cos(1000 * R) / R - 1e6 * np.sin(1000 * R),
    "sin(Z)*sin(R)": lambda R, Z: np.sin(Z) * np.cos(R) / R - 2 * np.sin(Z) * np.sin(R),
    "sin(Z)": lambda R, Z: -np.sin(Z),
    "sin(Z*10)": lambda R, Z: -np.sin(Z * 10) * 100,
    "sin(Z*100)": lambda R, Z: -np.sin(Z * 100) * 1e4,
    "sin(Z*1000)": lambda R, Z: -np.sin(Z * 1000) * 1e6,
    "Z": lambda R, Z: Z * 0,
    "r": r_func,
    "sin(theta)": sin_theta,
}
ana = dict()
ana["bracket(a, f)"] = {
    "R, Z": lambda R, Z: -1 / R,
    "R, R": lambda R, Z: 0 * R,
    "Z, R": lambda R, Z: 0 * R,
    "Z, Z": lambda R, Z: 0 * R,
    "sin(R), sin(Z)": lambda R, Z: -1 / R * np.cos(R) * np.cos(Z),
}
for a in 10, 100, 1000:
    ana["bracket(a, f)"][f"sin(R*{a}), sin(Z*{a})"] = (
        lambda R, Z, a=a: -1 / R * np.cos(a * R) * np.cos(a * Z) * a * a
    )
ana["bracket(a, f, OLD)"] = ana["bracket(a, f)"]
ana["FCI::Div_a_Grad_perp(a, f)"] = {
    "R, Z": lambda R, Z: 0 * R,
    "R, R": lambda R, Z: 0 * R + 2,
    "Z, R": lambda R, Z: Z / R,
    "Z, Z": lambda R, Z: 0 * Z + 1,
    "sin(R), sin(Z)": lambda R, Z: -np.sin(R) * np.sin(Z),
}
for a in 10, 100, 1000:
    ana["FCI::Div_a_Grad_perp(a, f)"][f"sin(R*{a}), sin(Z*{a})"] = (
        lambda R, Z, a=a: -np.sin(a * R) * np.sin(a * Z) * a * a
    )
ana["FCI::dagp(f)"] = ana["FCI::Div_a_Grad_perp(a, f)"]
ana["FCI::dagp_fv(f)"] = ana["FCI::Div_a_Grad_perp(a, f)"]


def get_ana(method, func):
    try:
        dic = ana[method]
    except KeyError:
        dic = ana_default
    try:
        return dic[func]
    except:
        print(method, func)
        raise


def doit(path):
    
    s = slice(2, -2), slice(None), slice(None)
    datasets = []

    # open the series of "BOUT.mesh_{m}.0.nc" files, 
    # saving them in a list `datasets`
    # and stop when no further files are found
    # save the counter variable `m`
    m = 0
    while 1:
        try:
            datasets += [ open_boutdataset(f'BOUT.mesh_{m}.0.nc', inputfilepath='BOUT.inp', geometry="toroidal",gridfilepath=f'circ_grid_{m}.nc',keep_yboundaries=False, is_mms_dump=True) ]
            m += 1
        except OSError:
            break
    # create a range from the number of files found
    nmesh = m
    meshrange = range(nmesh)

    # for each file, extract Z and R variables from
    # the previously saved datasets
    Zs = [collectvar(datasets, "Z", m) for m in meshrange]
    Rs = [collectvar(datasets, "R", m) for m in meshrange]
    
    # check that more than one dataset is supplied for the test
    assert len(Rs) > 1
    toplot = []

    # count the number of output variables in the BOUT.mesh_0.0.nc dataset
    # retain this variable `nvariable` to use in the loop over functions below
    nvariable = 0
    while 1:
        try:
            collectvar(datasets, f"out_{nvariable}")
        except KeyError:
            break
        nvariable += 1
    print(f"Checking {nvariable} variables for {nmesh} meshes in dir {path}")

    def defdiclist():
        return collections.defaultdict(list)
    nvariable = 1
    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for i in range(nvariable):
        l2norm = []
        nylist = []
        for m in meshrange:
            o = collectvar(datasets, f"out_{i}", m)
            attrs = o.attrs
            ops, inp = attrs["operator"], attrs["inp"]
            a = get_ana(ops, inp)(Rs[m], Zs[m])
            e = (o - a)[s]

            thisl2 = np.sqrt(np.mean(e**2))
            l2norm.append(thisl2)
            out[inp][ops].append(thisl2)
            nylist.append(Rs[m].shape[1])
        
        if not np.any(a):
            print(ops, inp)
            continue

        ord = []
        for i0 in range(len(l2norm) - 1):
            a, b = nylist[i0 : i0 + 2]
            dx = b / a
            a, b = l2norm[i0 : i0 + 2]
            de = a / b
            ord += [np.log(de) / np.log(dx)]
        if not np.isclose(
            ord[-1], get_convergence(attrs["operator"]), atol=0.25, rtol=0
        ):
            state = "❌"

            global success, failed, failed2
            failed[ops].append(inp)
            if ops not in expected_failures:
                success = False
                failed2.add(ops)
        else:
            state = "✅"
        print(
            state,
            i,
            np.array(l2norm),
            np.array(ord),
            {k: v for k, v in attrs.items() if "_" not in k and v},
        )

        toplot.append((attrs["inp"], attrs["operator"], nylist, l2norm))
        label = f'{attrs["inp"]} {attrs["operator"]}'
        with open(f"result_real_{i}.txt", "w") as f:
            f.write("real\n")
            f.write(f"{label}\n")
            f.write(" ".join([str(x) for x in nylist]))
            f.write("\n")
            f.write(" ".join([str(x) for x in l2norm]))
            f.write("\n")
    toplot2 = dict()
    for a, b, c, d in toplot:
        toplot2[a] = []
        # toplot2[b] = []
    for a, b, c, d in toplot:
        toplot2[a].append((b, c, d))
        # toplot2[b].append((a, c, d))
    out = {k: dict(v) for k, v in out.items()}
    #with open(f"{path}/l2_data.pkl", "wb") as f:
    #    pickle.dump(out, f, pickle.HIGHEST_PROTOCOL)
    if 1:
        for k, vs in toplot2.items():
            plt.figure()
            for ab, c, d in vs:
                plt.plot(c, d, "x-", label=ab)
            plt.title(k)
            plt.legend()
            plt.gca().set_yscale("log")
            plt.gca().set_xscale("log")
        plt.show()


if __name__ == "__main__":
    # to run in the present directory
    # $ python3 mms.py ./
    print(sys.argv)
    if sys.argv[1:]:
        args = sys.argv[1:]
    else:
        args = glob.iglob("mms*/")
    print(args)
    for p in args:
        print(p)
        doit(p)

    if failed:
        print(failed)
    if failed2:
        print(failed2)
    sys.exit(0 if success else 1)
