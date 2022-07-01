# Calculate cross-coherence between two signals

import numpy as np
import matplotlib.pyplot as plt

def cross_coherence(signal1, signal2):
    assert len(signal1) == len(signal2)
    
    # Normalise signals
    norm1 = ((signal1 - np.mean(signal1)) / np.std(signal1)).flatten()
    norm2 = ((signal2 - np.mean(signal2)) / np.std(signal2)).flatten()

    # Make a 2D histogram
    nbins = int(len(signal1)**(1./4) + 0.5)

    return np.histogram2d(norm1, norm2, bins=nbins, density=True, range=[[-2, 2], [-2, 2]])

def plot_cross_coherence(signal1, signal2, axis=None):
    """ Make a 2D cross-coherence plot of two signals
    
    The plot is a 2D histogram of the normalised signals
    i.e removing the mean, and dividing by the standard deviation

    signal1 will be the horizontal axis
    signal2 is the vertical axis

    axis is an optional axis to use. If not given then a
    new figure and axis is created.

    Returns the plot axis
    """
    hist, xedges, yedges = cross_coherence(signal1, signal2)

    if axis is None:
        fig, axis = plt.subplots()

    axis.set_aspect("equal")
    X, Y = np.meshgrid(xedges, yedges)
    axis.pcolormesh(X, Y, hist)

    Xc, Yc = np.meshgrid(0.5*(xedges[0:-1] + xedges[1:]), 0.5*(yedges[0:-1] + yedges[1:]))
    histmax = np.max(hist)
    axis.contour(Xc, Yc, hist, levels=[0.3*histmax, 0.5*histmax, 0.8*histmax])

    axis.axvline(0.0, linestyle='--', color='k')
    axis.axhline(0.0, linestyle='--', color='k')

    return axis

def acdc(f):
    dc = np.mean(f, axis=-1)
    ac = np.copy(f)
    for z in range(ac.shape[-1]):
        ac[...,z] -= dc
    return ac, dc

if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Directory containing the BOUT++ data")
    parser.add_argument("-s", "--fontsize", type=int, default=14, help="Font size, default 14")
    parser.add_argument('-x', "--nox", action='store_true', help="Don't use X display. Uses Agg backend.")
    parser.add_argument('-t0', "--first_time", type=int, default=0, help="First time index")
    parser.add_argument('-t1', "--last_time", type=int, default=-1, help="Last time index")
    parser.add_argument('-x0', "--first_x", type=int, default=0, help="First X index")
    parser.add_argument('-x1', "--last_x", type=int, default=-1, help="Last X index")
    args = parser.parse_args()
    
    import matplotlib
    matplotlib.rcParams.update({'font.size': args.fontsize})
    if args.nox:
        matplotlib.use('Agg')

    path = args.path

    tind = [args.first_time, args.last_time]
    xind = [args.first_x, args.last_x]
    
    # Read data
    from boutdata import collect
    n = collect("ne", path=path, tind=tind, xind=xind)
    phi = collect("phi", path=path, tind=tind, xind=xind)

    # Average Z direction to get profiles
    nac, ndc = acdc(n)
    phiac, phidc = acdc(phi)
    
    plot_cross_coherence(phiac.flatten(), nac.flatten())
    plt.xlabel(r"$A\left(\phi\right)/\sigma$")
    plt.ylabel(r"$A\left(n_e\right)/\sigma$")
    plt.title("Cross coherence")

    plt.savefig("cross_coherence.png")
    plt.savefig("cross_coherence.pdf")
    
    plt.show()

