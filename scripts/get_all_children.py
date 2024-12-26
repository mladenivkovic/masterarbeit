#!/usr/bin/env python3

#--------------------------------------------------
# Find and print out all subhaloes
# for a given halo
#
# Usage: get_all_children.py output_XXXXX clumpID
#--------------------------------------------------

import numpy as np
import os
from sys import argv


# ======================
class clumpdata():
# ======================
    """
    Data from clump_XXXXX.txtYYYYY
    """

    def __init__(self, par):
        """
        par: params object
        """

        # read-in
        self.clumpid = np.zeros(1)  # clump ID
        self.parent = np.zeros(1)  # parent ID
        self.level = np.zeros(1)  # clump level

        return

    # ----------------------------------
    def read_clumpdata(self, par):
    # ----------------------------------
        """
        Reads in the clump data.
        Only for the z = 0 directory.
        par: params object
        """

        raw_data = [None for i in range(par.ncpu)]
        dirnrstr = str(par.snapnr).zfill(5)
        dirname = par.srcdir

        i = 0
        for cpu in range(par.ncpu):
            fname = os.path.join(dirname, 'clump_' + dirnrstr + '.txt' + str(cpu + 1).zfill(5))
            new_data = np.loadtxt(fname, dtype='int', skiprows=1, usecols=[0, 1, 2])
            if new_data.ndim == 2:
                raw_data[i] = new_data
                i += 1
            elif new_data.shape[0] == 3:  # if only 1 row is present in file
                raw_data[i] = np.atleast_2d(new_data)
                i += 1

        fulldata = np.concatenate(raw_data[:i], axis=0)

        self.clumpid = fulldata[:, 0]
        self.level = fulldata[:, 1]
        self.parent = fulldata[:, 2]

        return


    # ----------------------------------
    def find_children(self, haloID):
    # ----------------------------------
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
            for i, cid in enumerate(self.clumpid):
                if self.parent[i] in this_level_parents and cid != haloID:
                    last_added.append(cid)

            if len(last_added) == 0:
                break

            if loopcounter == 100:
                print("Finished 100 iterations, we shouldn't be this deep")
                break

        return children[1:]



# ===================
class params():
# ===================
    """
    Global parameters to be stored
    """

    def __init__(self):
        self.workdir = ""  # current work directory
        self.srcdir = "" # output_XXXXX directory
        self.snapnr = 0 # XXXXX from output_XXXXX
        self.ncpu = 1
        self.clumpid = 0  # which clump ID to work for.

        return

    # -----------------------------
    def read_cmdlineargs(self):
    # -----------------------------
        """
        Reads in the command line arguments and stores them in the
        global_params object.
        """

        nargs = len(argv)
        i = 1  # first cmdlinearg is filename of this file, so skip it

        outputdir = argv[1]
        if not os.path.exists(outputdir):
            print("Didn't find directory", outputdir)
            quit()
        abspath = os.path.abspath(outputdir)
        prefix, outputXXXXX = os.path.split(abspath)

        self.srcdir = outputdir
        self.snapnr = int(outputXXXXX[-5:])

        self.clumpid = int(argv[2])

        return

    # --------------------------
    def get_output_info(self):
    # --------------------------
        """
        Read in the output info based on the files in the current
        working directory.
        Reads in last directory, ncpu, noutputs. Doesn't read infofiles.
        """

        self.workdir = os.getcwd()
        filelist = os.listdir(self.srcdir)

        outputlist = []
        for filename in filelist:
            if filename.startswith('clump_'):
                outputlist.append(filename)

        if len(outputlist) < 1:
            print("I didn't find any clump_XXXXX files in current working directory.")
            print("Are you in the correct workdir?")
            quit()

        outputlist.sort()

        # read ncpu from infofile in last output directory
        infofile = self.srcdir + '/' + 'info_' + str(self.snapnr).zfill(5) + '.txt'
        f = open(infofile, 'r')
        ncpuline = f.readline()
        line = ncpuline.split()

        self.ncpu = int(line[-1])

        f.close()

        return




# ===================================
if __name__ == '__main__':
# ===================================

    p = params()

    # Read cmdlineargs, available output, get global parameters
    p.read_cmdlineargs()
    p.get_output_info()

    cd = clumpdata(p)
    cd.read_clumpdata(p)

    children = cd.find_children(p.clumpid)

    children.sort()

    #  i = 0
    #  for c in children:
    #      print("{0:8d}, ".format(c), end="")
    #      i += 1
    #      if i == 10:
    #          print("&")
    #          i = 0
    #  print()
    #  print(len(children))

    for c in children:
        print("{0:12d}".format(c))

    print("Halo  {0:12d}".format(p.clumpid))
    print("Total {0:12d}".format(len(children)))
