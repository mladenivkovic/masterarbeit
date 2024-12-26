#!/usr/bin/env python3

#-----------------------------------------------
# Read in output of get_clump_ID_of_orphans,
# and create a list of main haloes that all the
# orphan galaxies are in.
#
# You first need to compile and run
# get_clump_ID_of_orphans!
#
# usage:
# get_orphan_halos.py <outputnr>
#-----------------------------------------------

import os
import sys
import numpy as np

# only work for orphans with mass above threshold below
mthresh = 1e10



def read_cmdilneargs():
    """
    Read in directory number,
    generate filename to read in.
    """

    try:
        dirnr = int(sys.argv[1])
    except ValueError :
        abort('Argument was not an integer')

    orphanfile = 'orphan_clump_IDs-'+str(dirnr).zfill(5)+".txt"
    if not os.path.exists(orphanfile):
        print("File", orphanfile, "not found.")
        quit()

    return dirnr, orphanfile

def read_orphan_data(orphanfile):
    """
    Read in clump IDs of orphans and their mass
    """

    clumpID = np.loadtxt(orphanfile, dtype=np.int, skiprows=1, usecols=[1])
    orph_tot = clumpID.shape[0]
    mass = np.loadtxt(orphanfile, dtype=np.float, skiprows=1, usecols=[5])

    # Filter orphan data
    mask = clumpID > 0
    clumpID = clumpID[mask]
    mass = mass[mask]

    mask = mass > mthresh
    clumpID = clumpID[mask]
    mass = mass[mask]

    print("keeping", clumpID.shape[0], "/", orph_tot, "orphans") 

    return clumpID, mass




def read_clump_data(dirnr):
    """
    Read in clump data 
    """

    dirnrstr = str(dirnr).zfill(5)
    srcdir = "output_"+dirnrstr
    workdir = os.getcwd()
    filelist = os.listdir(srcdir)

    # read ncpu from infofile in last output directory
    infofile = srcdir + '/' + 'info_' + dirnrstr + '.txt'
    f = open(infofile, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()
    ncpu = int(line[-1])
    f.close()


    # Now read in clump data file by file
    raw_int_data = [None for i in range(ncpu)]
    raw_float_data = [None for i in range(ncpu)]

    i = 0
    for cpu in range(ncpu):
        fname = os.path.join(srcdir, 'clump_' + dirnrstr + '.txt' + str(cpu + 1).zfill(5))
        new_int_data = np.loadtxt(fname, dtype='int', skiprows=1, usecols=[0, 1, 2])
        new_float_data = np.loadtxt(fname, dtype='float', skiprows=1, usecols=[4, 5, 6])
        if new_int_data.ndim == 2:
            raw_int_data[i] = new_int_data
            raw_float_data[i] = new_float_data
            i += 1
        elif new_int_data.shape[0] == 3:  # if only 1 row is present in file
            raw_int_data[i] = np.atleast_2d(new_int_data)
            raw_float_data[i] = np.atleast_2d(new_float_data)
            i += 1

    fullintdata = np.concatenate(raw_int_data[:i], axis=0)
    fullfloatdata = np.concatenate(raw_float_data[:i], axis=0)

    clumpid = fullintdata[:, 0]
    level = fullintdata[:, 1]
    parent = fullintdata[:, 2]
    coords = fullfloatdata

    return clumpid, level, parent, coords



def find_main_halos(orphanClumpID, clumpID, clumpLevel, clumpParent, clumpCoords):
    """
    Find all unique main halos of orphans.
    """

    orphanHalos = []
    orphanHaloCoords = []
    checked = [[] for i in range(clumpLevel.max()+1)]

    for ocid in orphanClumpID:
        cid = ocid
        ind = np.where(clumpID == cid)[0]
        if ind.shape[0] != 1:
            print("potential error for clump", cid)
            continue
        lev = clumpLevel[ind].item()
        if cid in checked[lev]:
            continue
        else:
            checked[lev].append(cid)

        is_halo = (clumpID[ind].item() == clumpParent[ind].item())

        emergency_break = 0
        add = True
        while not is_halo:
            cid = clumpParent[ind].item()
            ind = np.where(clumpID == cid)
            lev = clumpLevel[ind].item()

            if cid in checked[lev]:
                add = False
                break
            else:
                checked[lev].append(cid)

            is_halo = (cid == clumpParent[ind].item())
            emergency_break += 1
            if emergency_break > 20:
                print("Looped over 20 clump levels?")
                quit()

        if add:
            orphanHalos.append(cid)
            orphanHaloCoords.append(clumpCoords[ind][0])

    return orphanHalos, orphanHaloCoords



def write_results(halos, coords, dirnr):
    """
    write results to file
    """

    fname = "unique_orphan_main_halos-"+str(dirnr).zfill(5)+".txt"
    f = open(fname, 'w')
    for i in range(len(halos)):
        line = "{0:16d} {1:12.8f} {2:12.8f} {3:12.8f}\n".format(
                    halos[i], coords[i][0], coords[i][1], coords[i][2]
                )
        f.write(line)

    return
 



def main():
    dirnr, orphanfile = read_cmdilneargs()
    orphanClumpID, orphanMass = read_orphan_data(orphanfile)
    clumpID, clumpLevel, clumpParent, clumpCoords = read_clump_data(dirnr)
    orphanHalos, orphanHaloCoords = find_main_halos(orphanClumpID, clumpID, clumpLevel, clumpParent, clumpCoords)

    print("keeping", len(orphanHalos), "main halos")

    if np.unique(orphanHalos).shape[0] != len(orphanHalos):
        print("You have non-unique main halos in there?")

    write_results(orphanHalos, orphanHaloCoords, dirnr)



if __name__ == "__main__" :
    main()




def find_children(haloID, clumpID, clumpParent):
    """
    Find the children for given clump ID.
    haloID: clump ID for which to work for

    returns:
        children:   list of children IDs of clumpid
    """

    import copy

    children = []
    last_added = [haloID]

    loopcounter = 0
    while True:
        loopcounter += 1
        this_level_parents = copy.copy(last_added)
        children += this_level_parents
        last_added = []
        for i, cid in enumerate(clumpID):
            if clumpParent[i] in this_level_parents and cid != haloID:
                last_added.append(cid)

        if len(last_added) == 0:
            break

        if loopcounter == 100:
            print("Finished 100 iterations, we shouldn't be this deep")
            break


    return children[1:]


