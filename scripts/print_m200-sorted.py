#!/usr/bin/env python3


#-----------------------------------------------
# Look for all radial-profile-metadata-*dat
# files in given output directory, and print
# halo IDs sorted by M200c mass
#
# usage:
#   print_m200-sorted.py <outputnr>
# the script expects the directory output_<outputnr>
# to exists in the current workdir
#-----------------------------------------------

import sys
import os
import numpy as np

#  print_additional_details = True
print_additional_details = False
use_partcenter = True
print_only = 10


def get_output_dir():
    """
    Read in and construct the directory name to work with
    """
    
    try:
        outputnr = sys.argv[1]
    except IndexError:
        print("You need to give me the output directory number as a cmdline arg")
        quit(1)

    dirname = "output_"+str(outputnr).zfill(5)
    if not os.path.exists(dirname):
        print("Couldn't find directory", dirname)
        quit()

    return dirname

    
def get_filelist(outputdir):
    """
    Get all files that we want to read in
    """

    fulllist = os.listdir(outputdir)
    filelist = []
    if use_partcenter:
        for f in fulllist:
            if f.startswith("radial-profile-partcenter-metadata"):
                filelist.append(f)
    else:
        for f in fulllist:
            if f.startswith("radial-profile-metadata"):
                filelist.append(f)

    if len(filelist) == 0:
        print("Found no radial-profile-metadata* files")
        quit()

    return filelist



def read_file(filename):
    """
    Read the actual file
    """

    split = filename.split("-")
    halo = split[-1]
    halo = halo[:-4] # remove .dat
    halo = int(halo)

    with open(filename) as handle:
        for i in range(5):
            handle.readline()
        radius = float(handle.readline())
        for i in range(3):
            handle.readline()
        mass = float(handle.readline())

    return halo, mass, radius

 
def read_partcenter_file(filename):
    """
    Read the actual file
    """

    split = filename.split("-")
    halo = split[-1]
    halo = halo[:-4] # remove .dat
    halo = int(halo)

    with open(filename) as handle:
        for i in range(5):
            handle.readline()
        radius = float(handle.readline())
        for i in range(3):
            handle.readline()
        mass = float(handle.readline())
        for i in range(3):
            handle.readline()
        unit_l_Mpc = float(handle.readline())
        radius *= unit_l_Mpc

    return halo, mass, radius

   


def read_data(outputdir, filelist):
    """
    Read in the actual data
    """

    halo = []
    mass = []
    radius = []

    for f in filelist:
        fname = os.path.join(outputdir, f)
        if use_partcenter:
            h, m, r = read_partcenter_file(fname)
        else:
            h, m, r = read_file(fname)
        halo.append(h)
        mass.append(m)
        radius.append(r)

    halo = np.array(halo)
    mass = np.array(mass)
    radius = np.array(radius)

    sortind = np.argsort(mass)
    sortind = np.flip(sortind)
    halo = halo[sortind]
    mass = mass[sortind]
    radius = radius[sortind]

    return halo, mass, radius


def print_results(halo, mass, radius):

    if print_additional_details:
        for i, h in enumerate(halo):
            print("{0:12d} {1:12.3e} Msol {2:12.3e} Mpc".format(h, mass[i], radius[i]))
    else:
        counter = 0
        for h in halo:
            print("{0:15d}".format(h))
            counter += 1
            if print_only is not None and counter == print_only:
                break
    


def main():
    outputdir = get_output_dir()
    filelist = get_filelist(outputdir)
    halo, mass, radius = read_data(outputdir, filelist)
    print_results(halo, mass, radius)

if __name__ == "__main__":
    main()
