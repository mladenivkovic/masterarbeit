#!/usr/bin/env python3

#--------------------------------------------------
# Trace particles back in time over multiple
# snapshots. Use extract_clumpparticles.py first
# to extract the particles and boundary.
#
# usage: 
#   plot_clumpparticles.py clumpparticles_XXXX.pkl
# where clumpparticles_*pkl is the output of 
# extract_clumpparticles.py
#--------------------------------------------------


nsnapshots = 10 # how many snapshots to go back in time


import sys, os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from scipy.io import FortranFile
import pickle

Mpc    = 3.0857e+24 # Mpc in cm


def get_cmdlineargs():
    """
    Read in cmdline arguments.
    """

    try:
        clumppartfile = sys.argv[1]
    except:
        print("Usage: plot_clumpparticles.py clumpparticles_XXXXX-clump-YYY.pkl")
        quit()

    if not os.path.exists(clumppartfile):
        print("Didn't find file", clumppartfile)

    return clumppartfile


def get_particles(outputnr, clumppartIDs, ncpu, boundaries):
    """
    Read in particles data
    """

    print("reading")

    dirname = "output_{0:05d}".format(outputnr)
    partfile_base = "part_{0:05d}.out".format(outputnr)
    partfile_base = os.path.join(dirname, partfile_base)

    if not os.path.exists(partfile_base+"00001"):
        print("Couldn't find "+partfile_base+"00001 in", dirname)
        print("To plot particles, I require the particle output.")
        quit()


    all_particles = []
    clumpparticles = []
    nonclumpparticles = []

    xmin, xmax, ymin, ymax, zmin, zmax = boundaries

    for cpu in range(ncpu):
        partfile = partfile_base+str(cpu+1).zfill(5)
        pf = FortranFile(partfile, "r")

        _ = pf.read_ints() # ncpu
        _ = pf.read_ints() # ndim
        _ = pf.read_ints() # nparts
        _ = pf.read_ints() # localseed
        _ = pf.read_ints() # nstar_tot
        _ = pf.read_reals('d') # mstar_tot
        _ = pf.read_reals('d') # mstar_lost
        _ = pf.read_ints() # nsink


        x = pf.read_reals('d')
        y = pf.read_reals('d')
        z = pf.read_reals('d')

        _ = pf.read_reals('d') # vx
        _ = pf.read_reals('d') # vz
        _ = pf.read_reals('d') # vz
        _ = pf.read_reals('d') # mass
        idp = pf.read_ints()

        # add periodicity corrections
        if xmin < 0:
            x[x>0.5] -= 1.
        if xmax > 1.:
            x[x<0.5] += 1.
        if ymin < 0:
            y[y>0.5] -= 1.
        if ymax > 1.:
            y[y<0.5] += 1.
        if zmin < 0:
            z[z>0.5] -= 1.
        if zmax > 1.:
            z[z<0.5] += 1.

        for i in range(x.shape[0]):
            in_box = False
            is_clumpparticle = False

            in_x = x[i] >= xmin and x[i] <= xmax
            in_y = y[i] >= ymin and y[i] <= ymax
            in_z = z[i] >= zmin and z[i] <= zmax
            if in_x and in_y and in_z:
                all_particles.append([x[i], y[i], z[i]])
                in_box = True


            if idp[i] in clumppartIDs:
                clumpparticles.append([x[i], y[i], z[i]])
                is_clumpparticle = True

                if x[i] < xmin:
                    print("Found clumpparticle outside radius (1), resize", x[i], xmin)
                if x[i] > xmax:
                    print("Found clumpparticle outside radius (2), resize", x[i], xmax)
                if y[i] < ymin:
                    print("Found clumpparticle outside radius (3), resize", y[i], ymin)
                if y[i] > ymax:
                    print("Found clumpparticle outside radius (4), resize", y[i], ymax)
                if z[i] < zmin:
                    print("Found clumpparticle outside radius (5), resize", z[i], zmin)
                if z[i] > zmax:
                    print("Found clumpparticle outside radius (6), resize", z[i], zmax)

            if in_box and not is_clumpparticle:
                nonclumpparticles.append([x[i], y[i], z[i]])

        pf.close()

    all_particles = np.array(all_particles)
    clumpparticles = np.array(clumpparticles)
    nonclumpparticles = np.array(nonclumpparticles)


    return all_particles, clumpparticles, nonclumpparticles


def plot_particles(outputnr, outputnr_start, clumpID, allparts, clumpparts, nonclumpparts, boundaries):
    """
    Create the plot.
    """
    print("plotting")

    xmin, xmax, ymin, ymax, zmin, zmax = boundaries
    xrange = [xmin, xmax]
    yrange = [ymin, ymax]
    zrange = [zmin, zmax]
    nbins = 400

    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(3, 3, 1)
    ax2 = fig.add_subplot(3, 3, 2)
    ax3 = fig.add_subplot(3, 3, 3)

    ax4 = fig.add_subplot(3, 3, 4)
    ax5 = fig.add_subplot(3, 3, 5)
    ax6 = fig.add_subplot(3, 3, 6)

    ax7 = fig.add_subplot(3, 3, 7)
    ax8 = fig.add_subplot(3, 3, 8)
    ax9 = fig.add_subplot(3, 3, 9)

    hist1, _, _ = np.histogram2d(nonclumpparts[:,0], nonclumpparts[:,1], range=(xrange, yrange), bins=nbins)
    hist1 += 1
    hist2, _, _ = np.histogram2d(nonclumpparts[:,1], nonclumpparts[:,2], range=(yrange, zrange), bins=nbins)
    hist2 += 1
    hist3, _, _ = np.histogram2d(nonclumpparts[:,0], nonclumpparts[:,2], range=(xrange, zrange), bins=nbins)
    hist3 += 1

    hist7, _, _ = np.histogram2d(allparts[:,0], allparts[:,1], range=(xrange, yrange), bins=nbins)
    hist7 += 1
    hist8, _, _ = np.histogram2d(allparts[:,1], allparts[:,2], range=(yrange, zrange), bins=nbins)
    hist8 += 1
    hist9, _, _ = np.histogram2d(allparts[:,0], allparts[:,2], range=(xrange, zrange), bins=nbins)
    hist9 += 1

    imshow_kwargs = {
        "origin": 'lower', 
        "norm": colors.LogNorm(),
        "cmap": "Blues"
    }

    # Plot background
    ax1.imshow(hist1.T, extent=(xmin, xmax, ymin, ymax), **imshow_kwargs)
    ax2.imshow(hist2.T, extent=(ymin, ymax, zmin, zmax), **imshow_kwargs)
    ax3.imshow(hist3.T, extent=(xmin, xmax, zmin, zmax), **imshow_kwargs)

    ax4.imshow(hist1.T, extent=(xmin, xmax, ymin, ymax), **imshow_kwargs)
    ax5.imshow(hist2.T, extent=(ymin, ymax, zmin, zmax), **imshow_kwargs)
    ax6.imshow(hist3.T, extent=(xmin, xmax, zmin, zmax), **imshow_kwargs)

    # Plot background
    ax7.imshow(hist7.T, extent=(xmin, xmax, ymin, ymax), **imshow_kwargs)
    ax8.imshow(hist8.T, extent=(ymin, ymax, zmin, zmax), **imshow_kwargs)
    ax9.imshow(hist9.T, extent=(xmin, xmax, zmin, zmax), **imshow_kwargs)

    scatter_kwargs = {
        "alpha": 0.5,
        "lw": 0,
        "c": "red", 
        "s": 1, 
        "marker": "."
    }

    ax1.scatter(clumpparts[:,0], clumpparts[:,1], **scatter_kwargs)
    ax2.scatter(clumpparts[:,1], clumpparts[:,2], **scatter_kwargs)
    ax3.scatter(clumpparts[:,0], clumpparts[:,2], **scatter_kwargs)


    for ax in [ax1, ax4, ax7]:
        ax.set_xlim(xrange)
        ax.set_ylim(yrange)
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    for ax in [ax2, ax5, ax8]:
        ax.set_xlim(yrange)
        ax.set_ylim(zrange)
        ax.set_xlabel("y")
        ax.set_ylabel("z")
    ax2.set_title("Clumpparticles+Background")
    ax5.set_title("Non-Clumpparticles")
    ax8.set_title("All particles")

    for ax in [ax3, ax6, ax9]:
        ax.set_xlim(xrange)
        ax.set_ylim(zrange)
        ax.set_xlabel("x")
        ax.set_ylabel("z")

    fig.suptitle("output_{0:05d} clump {1:d}@output_{2:05d}".format(outputnr, clumpID, outputnr_start))
    fig.tight_layout()

    #  plt.show()
    plt.savefig("clump_particle_tracing_output_{0:05d}-{1:d}.png".format(outputnr, clumpID), dpi=240)


def read_dump(clumppartfile):
    """
    Read pickle file
    """
    
    pf = open(clumppartfile, "rb")
    outputnr = pickle.load(pf)
    clumpID = pickle.load(pf)
    ncpu = pickle.load(pf)
    boundaries = pickle.load(pf)
    clumpparticles = pickle.load(pf)
    pf.close()

    return outputnr, clumpID, ncpu, boundaries, clumpparticles


if __name__ == "__main__":

    clumppartfile = get_cmdlineargs()

    outputnr, clumpID, ncpu, boundaries, clumpparticles = read_dump(clumppartfile)

    for snap in range(nsnapshots):
        output = outputnr - snap
        allparts, clumpparts, nonclumpparts = get_particles(output, clumpparticles, ncpu, boundaries)
        plot_particles(output, outputnr, clumpID, allparts, clumpparts, nonclumpparts, boundaries)

    print("Done")

