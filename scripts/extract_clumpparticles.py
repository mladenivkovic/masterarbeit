#!/usr/bin/env python3

#--------------------------------------------------
# Extract particle IDs of a given clump
# and dump them into a pickle file.
# The pickle file is intended to be read in
# by plot_traced_clumpparticles.py
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
import pickle

Mpc    = 3.0857e+24 # Mpc in cm

def get_cmdlineargs():
    """
    Read in cmdline arguments.
    """

    try:
        outputnr = int(sys.argv[1])
        clumpID = int(sys.argv[2])
    except:
        print("Usage: extract_clumpparticles.py outputnr clump-ID")
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


def get_clumpparticles(outputnr, clumpID, ncpu, boundaries):
    """
    Read in particles data, and dump the particle IDs
    """

    dirname = "output_{0:05d}".format(outputnr)
    partfile_base = "part_{0:05d}.out".format(outputnr)
    partfile_base = os.path.join(dirname, partfile_base)
    unbindfile_base = "unbinding_{0:05d}.out".format(outputnr)
    unbindfile_base = os.path.join(dirname, unbindfile_base)

    if not os.path.exists(unbindfile_base+"00001"):
        print("Couldn't find "+unbindfile_base+"00001 in", dirname)
        print("To plot particles, I require the unbinding output.")
        quit()


    clumpparticles = []

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
        _ = pf.read_reals('d') # x
        _ = pf.read_reals('d') # y
        _ = pf.read_reals('d') # z
        _ = pf.read_reals('d') # vx
        _ = pf.read_reals('d') # vy
        _ = pf.read_reals('d') # vz
        _ = pf.read_reals('d') # mass

        idp = pf.read_ints() # IDs

        unbfile = unbindfile_base+str(cpu+1).zfill(5)
        unbffile = FortranFile(unbfile, "r")
        clumpidp = unbffile.read_ints()

        for i in range(clumpidp.shape[0]):

            if clumpID == abs(clumpidp[i]):
                clumpparticles.append(idp[i])

        pf.close()
        unbffile.close()

    clumpparticles = np.array(clumpparticles)

    print("Found", clumpparticles.shape[0], "particles")

    return clumpparticles


def dump_data(clumpparticles, boundaries, ncpu, outputnr, clumpID):
    """
    Dump pickle file
    """

    picklefilename = "clumpparticles_{0:05d}-clump-{1:d}.pkl".format(outputnr, clumpID)
    pf = open(picklefilename, "wb")

    pickle.dump(outputnr, pf)
    pickle.dump(clumpID, pf)
    pickle.dump(ncpu, pf)
    pickle.dump(boundaries, pf)
    pickle.dump(clumpparticles, pf)
    pf.close()





if __name__ == "__main__":

    outputnr, clumpID = get_cmdlineargs()
    ncpu, unit_l = read_infofile(outputnr)
    boundaries = get_clumpdata(outputnr, clumpID, ncpu, unit_l)
    clumpparticles = get_clumpparticles(outputnr, clumpID, ncpu, boundaries)
    dump_data(clumpparticles, boundaries, ncpu, outputnr, clumpID)
