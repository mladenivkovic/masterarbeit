#!/usr/bin/env python3

#---------------------------------------------
# Read in all galaxies from 
# output_XXXXX/galaxies_XXXXX.txtYYYYY files
# and histogram their masses for "regular"
# galaxies and for orphans
#
# Usage:
#   galaxy_masses_check_plot.py output_XXXXX
# or
#   galaxy_masses_check_plot.py XXXXX
#---------------------------------------------

import sys, os
import numpy as np
from matplotlib import pyplot as plt

ourputnrstr = ""


def get_output_info(srcdir = "."):
    """
    Read in the output info based on the files in the current
    working directory.
    Reads in last directory, ncpu, noutputs, workdir.

    returns:        ncpu
    """

    filelist = os.listdir(srcdir)
    infofile = None
    for f in filelist:
        if f.startswith("info_"):
            infofile = f
            break

    if infofile is None:
        print("Infofile is none?")
        quit()

    infofilepath = os.path.join(srcdir, infofile)
    f = open(infofilepath, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()

    ncpu = int(line[-1])
    f.close()

    print("got ncpu", ncpu)

    return ncpu


def get_srcdir():
    try:
        dirarg = sys.argv[1]
    except IndexError:
        print("I need either an output_XXXXX directory as a cmdline arg, or XXXXX as an integer")
        quit()
    try:
        dirint = int(dirarg)
        srcdir = "output_"+str(dirint).zfill(5)
    except ValueError:
        srcdir = dirarg

    if not os.path.exists(srcdir):
        print("Can't find srcdir=", srcdir)
        quit()

    print("Got srcdir", srcdir)

    global outputnrstr
    if srcdir[-1] == '/':
        outputnrstr = srcdir[-6:-1]
    else:
        outputnrstr = srcdir[-5:]

    return srcdir



def read_galaxy_masses(srcdir, ncpu):
    """
    Read in galaxy and orphan masses

    """

    filelist = os.listdir(srcdir)
    galaxyfiles = []
    for f in filelist:
        if f.startswith("galaxies_"):
            galaxyfiles.append(f)

    if len(galaxyfiles) != ncpu:
        print("number of galaxy files (", len(galaxyfiles), ") != ncpu: (", ncpu, ")")


    galaxy_masses = []
    orphan_masses = []
    for gf in galaxyfiles:
        galaxyfilepath = os.path.join(srcdir, gf)
        clump = np.loadtxt(galaxyfilepath, skiprows=1, usecols=[0], dtype=np.int)
        mass = np.loadtxt(galaxyfilepath, skiprows=1, usecols=[1], dtype=np.float)
        galaxy_masses.append(np.copy(mass[clump>0]))
        orphan_masses.append(np.copy(mass[clump==0]))

    gm = np.hstack(galaxy_masses)
    om = np.hstack(orphan_masses)

    return gm, om


def plot_masses(galm, orphm):
    """
    Create histograms of galaxy and orphan masses and plot them
    """

    print("Galaxy min/max:{0:.3E} {1:.3E}".format(galm.min(), galm.max()))
    print("Orphan min/max:{0:.3E} {1:.3E}".format(orphm.min(), orphm.max()))


    loggm = np.log10(galm)
    logom = np.log10(orphm)

    nbins = 200

    plt.figure()
    plt.hist(loggm, bins=nbins, label="galaxies with host", histtype="step", bottom=1, log=True)
    plt.hist(logom, bins=nbins, label="orphans", histtype="step", bottom=1, log=True)
    plt.xlabel(r"log $M_*/M_\odot$")
    plt.ylabel(r"$N + 1$")
    plt.legend()
    plt.grid()
    #  plt.show()
    plt.savefig("galaxy_masses_check-"+outputnrstr+".png", dpi=300)

    nt1 = np.count_nonzero(logom>9.)
    print("orphans > 10^9 M_Sol:", nt1)
    nt2 = np.count_nonzero(logom>10.)
    print("orphans > 10^10 M_Sol:", nt2)

    
 


def main():

    # get directory to work with
    srcdir = get_srcdir()
    ncpu = get_output_info(srcdir)

    gm, om = read_galaxy_masses(srcdir, ncpu)

    plot_masses(gm, om)

    



if __name__ == "__main__":
    main()
