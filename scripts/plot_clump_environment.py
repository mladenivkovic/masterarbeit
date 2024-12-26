#!/usr/bin/env python3


#--------------------------------------------------
# Plot a clump and all other particles
# in a given environment
#
# usage:
#   plot_clump_environment.py outputnr clump-ID
#--------------------------------------------------


radius = 2.5 # in  Mpc. Defines environment


import sys, os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from scipy.io import FortranFile

Mpc    = 3.0857e+24 # Mpc in cm


def get_cmdlineargs():
    """
    Read in cmdline arguments.
    """

    try:
        outputnr = int(sys.argv[1])
        clumpID = int(sys.argv[2])
    except:
        print("Usage: plot_clump_environment.py outputnr clump-ID")
        quit()

    outputdir = "output_{0:05d}".format(outputnr)
    if not os.path.exists(outputdir):
        print("Didn't find directory", outputdir)

    return outputnr, clumpID


def read_infofile(outputnr):
    """
    Read in data from the info file.
    outputnr: XXXXX in output_XXXXX
    """

    outputdir = "output_{0:05d}".format(outputnr)
    infofilename = "info_{0:05d}.txt".format(outputnr)


    infofile = os.path.join(outputdir, infofilename)

    f = open(infofile, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()
    ncpu = int(line[-1])

    # skip 14 lines
    for i in range(14):
        f.readline()

    unitlline = f.readline()
    line = unitlline.split()
    unit_l = float(line[-1])

    f.close()

    return ncpu, unit_l


def get_clumpdata(outputnr, clumpID, ncpu, unit_l):
    """
    Read in clump position and find boundaries
    based on "radius"
    """

    dirname = "output_{0:05d}".format(outputnr)
    clumpfile_base = "clump_{0:05d}.txt".format(outputnr)
    clumpfile_base = os.path.join(dirname, clumpfile_base)

    xpeak, ypeak, zpeak = None, None, None

    for cpu in range(ncpu):
        fname = clumpfile_base + str(cpu + 1).zfill(5)
        new_data = np.loadtxt(fname, dtype='float', skiprows=1, usecols=[0, 4, 5, 6])

        for c in range(new_data.shape[0]):
            if clumpID == new_data[c, 0]:
                xpeak, ypeak, zpeak = new_data[c, 1:4]
                # quit when you got what you need
                break
        
        # quit when you got what you need
        if xpeak is not None:
            break

    
    if xpeak is None:
        print("Didn't find data for clump", clumpID)
        quit()

    dx = radius * Mpc / unit_l
    xmin = xpeak - dx
    xmax = xpeak + dx
    ymin = ypeak - dx
    ymax = ypeak + dx
    zmin = zpeak - dx
    zmax = zpeak + dx


    return [xmin, xmax, ymin, ymax, zmin, zmax]


def get_particles(outputnr, clumpID, ncpu, boundaries):
    """
    Read in particles data
    """

    print("reading")

    dirname = "output_{0:05d}".format(outputnr)
    partfile_base = "part_{0:05d}.out".format(outputnr)
    partfile_base = os.path.join(dirname, partfile_base)
    unbindfile_base = "unbinding_{0:05d}.out".format(outputnr)
    unbindfile_base = os.path.join(dirname, unbindfile_base)

    if not os.path.exists(unbindfile_base+"00001"):
        print("Couldn't find "+unbindfile_base+"00001 in", dirname)
        print("To plot particles, I require the unbinding output.")
        quit()


    particles = []
    clumpparticles = []
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

        unbfile = unbindfile_base+str(cpu+1).zfill(5)
        unbffile = FortranFile(unbfile, "r")
        clumpidp = unbffile.read_ints()

        for i in range(clumpidp.shape[0]):
            in_x = x[i] >= xmin and x[i] <= xmax
            in_y = y[i] >= ymin and y[i] <= ymax
            in_z = z[i] >= zmin and z[i] <= zmax
            if in_x and in_y and in_z:
                particles.append([x[i], y[i], z[i]])

            if clumpID == abs(clumpidp[i]):
                clumpparticles.append([x[i], y[i], z[i]])

        pf.close()
        unbffile.close()

    particles = np.array(particles)
    clumpparticles = np.array(clumpparticles)

    return particles.astype(np.float32), clumpparticles.astype(np.float32)


def plot_particles(outputnr, clumpID, parts, clumpparts, boundaries):
    """
    Create the plot.
    """
    print("plotting")

    xmin, xmax, ymin, ymax, zmin, zmax = boundaries
    xrange = [xmin, xmax]
    yrange = [ymin, ymax]
    zrange = [zmin, zmax]
    nbins = 400

    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2)
    ax3 = fig.add_subplot(1, 3, 3)

    hist1, _, _ = np.histogram2d(parts[:,0], parts[:,1], range=(xrange, yrange), bins=nbins)
    hist1 += 1
    hist2, _, _ = np.histogram2d(parts[:,1], parts[:,2], range=(yrange, zrange), bins=nbins)
    hist2 += 1
    hist3, _, _ = np.histogram2d(parts[:,0], parts[:,2], range=(xrange, zrange), bins=nbins)
    hist3 += 1

    imshow_kwargs = {
        "origin": 'lower', 
        "norm": colors.LogNorm(),
        "cmap": "Blues"
    }

    ax1.imshow(hist1.T, extent=(xmin, xmax, ymin, ymax), **imshow_kwargs)
    ax2.imshow(hist2.T, extent=(ymin, ymax, zmin, zmax), **imshow_kwargs)
    ax3.imshow(hist3.T, extent=(xmin, xmax, zmin, zmax), **imshow_kwargs)

    scatter_kwargs = {
        "alpha": 0.6,
        "lw": 0,
        "c": "red", 
        "s": 1, 
        "marker": "."
    }

    ax1.scatter(clumpparts[:,0], clumpparts[:,1], **scatter_kwargs)
    ax2.scatter(clumpparts[:,1], clumpparts[:,2], **scatter_kwargs)
    ax3.scatter(clumpparts[:,0], clumpparts[:,2], **scatter_kwargs)


    ax1.set_xlim(xrange)
    ax1.set_ylim(yrange)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax2.set_xlim(yrange)
    ax2.set_ylim(zrange)
    ax2.set_xlabel("y")
    ax2.set_ylabel("z")
    ax3.set_xlim(xrange)
    ax3.set_ylim(zrange)
    ax3.set_xlabel("x")
    ax3.set_ylabel("z")
    fig.suptitle("output_{0:05d} clump {1:d}".format(outputnr, clumpID))
    fig.tight_layout()

    #  plt.show()
    plt.savefig("clump_environment_plot_output_{0:05d}-{1:d}.png".format(outputnr, clumpID), dpi=240)



if __name__ == "__main__":

    outputnr, clumpID = get_cmdlineargs()
    ncpu, unit_l = read_infofile(outputnr)
    boundaries = get_clumpdata(outputnr, clumpID, ncpu, unit_l)
    parts, clumpparts = get_particles(outputnr, clumpID, ncpu, boundaries)
    plot_particles(outputnr, clumpID, parts, clumpparts, boundaries)

