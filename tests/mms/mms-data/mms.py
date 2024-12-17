from boutdata import collect as boutcollect
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


def collect(path, var, mesh=0):
    print(var, path, mesh)
    return boutcollect(
        var, path=path, prefix=f"BOUT.mesh_{mesh}", strict=True, info=False
    )

def collectgrid(path, var, mesh=0):
    print(var, path, mesh)
    return boutcollect(
        var, path=path, prefix=f"guards_grid_{mesh}", strict=True, info=False
    )

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
    if 0:
        r = os.system(f"build-cont-opt/hermes-mms -d {path} -q -q -q")
        if r:
            os.system(f"build-cont-opt/hermes-mms -d {path}")
            raise RuntimeError("bout++ failed")

    s = slice(2, -2), slice(None), slice(None)

    Rs = []
    m = 0
    while 1:
        try:
            Rs += [collectgrid(path, "Rxy", m)]
            m += 1
        except OSError:
            break
    mids = range(m)
    Zs = [collectgrid(path, "Zxy", m) for m in mids]
    assert len(Rs) == len(Zs)
    assert len(Rs) > 1
    toplot = []
    mmax = 0
    while 1:
        try:
            collect(path, f"out_{mmax}")
        except ValueError:
            break
        mmax += 1
    print(f"Checking {mmax} variables for {len(Zs)} meshes in dir {path}")

    def defdiclist():
        return collections.defaultdict(list)

    out = collections.defaultdict(lambda: collections.defaultdict(list))
    for i in list(range(mmax)):
        l2 = []
        lst = []
        for m in mids:
            o = collect(path, f"out_{i}", m)
            attrs = o.attributes
            ops, inp = attrs["operator"], attrs["inp"]
            a = get_ana(ops, inp)(Rs[m], Zs[m])
            e = (o - a)[s]
            if "dagp_fv" in ops and 0:
                f, axs = plt.subplots(1, 3)
                s2 = slice(2, -2, None), 0, slice(None)
                print([x.shape for x in [Rs[m][s2], Zs[m][s2], o[s2]]])
                ax = axs[0]
                p = ax.pcolormesh(Rs[m][s2], Zs[m][s2], o[s2])
                ax.set_title(f"{inp} {ops} output[s2]")
                plt.colorbar(p, ax=ax)
                ax = axs[1]
                p = ax.pcolormesh(Rs[m][s2], Zs[m][s2], a[s2])
                ax.set_title(f"{inp} {ops} analytic[s2]")
                plt.colorbar(p, ax=ax)
                ax = axs[2]
                emax = np.max(np.abs((o - a)[s2]))
                p = ax.pcolormesh(
                    Rs[m][s2],
                    Zs[m][s2],
                    (o - a)[s2],
                    cmap=plt.get_cmap("bwr"),
                    vmax=emax,
                    vmin=-emax,
                )
                ax.set_title(f"{inp} {ops} error[s2]")
                plt.colorbar(p, ax=ax)

            thisl2 = np.sqrt(np.mean(e**2))
            l2.append(thisl2)
            out[inp][ops].append(thisl2)
            lst.append(Rs[m].shape[2])
        if not np.any(a):
            print(ops, inp)
            continue
        plt.show()

        ord = []
        for i0 in range(len(l2) - 1):
            a, b = lst[i0 : i0 + 2]
            dx = b / a
            a, b = l2[i0 : i0 + 2]
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
            np.array(l2),
            np.array(ord),
            {k: v for k, v in attrs.items() if "_" not in k and v},
        )

        toplot.append((attrs["inp"], attrs["operator"], lst, l2))
        label = f'{attrs["inp"]} {attrs["operator"]}'
        with open(f"result_real_{i}.txt", "w") as f:
            f.write("real\n")
            f.write(f"{label}\n")
            f.write(" ".join([str(x) for x in lst]))
            f.write("\n")
            f.write(" ".join([str(x) for x in l2]))
            f.write("\n")
    toplot2 = dict()
    for a, b, c, d in toplot:
        toplot2[a] = []
        # toplot2[b] = []
    for a, b, c, d in toplot:
        toplot2[a].append((b, c, d))
        # toplot2[b].append((a, c, d))
    out = {k: dict(v) for k, v in out.items()}
    with open(f"{path}/l2_data.pkl", "wb") as f:
        pickle.dump(out, f, pickle.HIGHEST_PROTOCOL)
    if 0:
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
