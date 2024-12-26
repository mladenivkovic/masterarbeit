#!/usr/bin/env python3

#-----------------------------------------------
# This script finds and plots the most massive
# halos in the last output directory
#
# optional cmd line arg: how many halos to print
# in descending mass order
#-----------------------------------------------



from os import listdir, path
import numpy as np
import argparse


def parse_args():
    """
    Parse cmdline args.
    """

    parser = argparse.ArgumentParser(
        description="""
        Print halo IDs of most massive haloes in directory.
        By default, the script finds the directory with redshift
        closest to zero.
        """
    )

    parser.add_argument("-n", "--nhalos", dest="nhalos", type=int, default=10, help="number of halos to look for")
    parser.add_argument("-d", "--dir", dest="srcdir", type=str, default="None", help="snapshot directory to look in")

    args = parser.parse_args()
    return  args.nhalos, args.srcdir


#============================
def read_redshifts():
#============================

    from os import getcwd
    from os import listdir

    workdir = getcwd()
    filelist = listdir(workdir)

    outputlist = []
    outdirtemp = 'output_'
    for filename in filelist:
        if filename.startswith(outdirtemp):
            outputlist.append(filename)

    if len(outputlist)<1:
        print("I didn't find any output_XXXXX directories in current working directory.")
        print("Are you in the correct workdir?")
        quit()

    outputlist.sort()
    noutput = len(outputlist)

    outputnrs = np.array(range(noutput))
    a_exp = np.zeros(noutput)

    for i,out in enumerate(outputlist):
        fn = out+'/info_'+out[-5:]+'.txt'
        try:
            infofile = open(fn)
            for j in range(9):
                infofile.readline() # skip first 9 lines

            # get expansion factor
            aline = infofile.readline()
            astring, equal, aval = aline.partition("=")
            afloat = float(aval)
            a_exp[i] = afloat

        except IOError: # If file doesn't exist
            print("Didn't find any info data in ", out)
            break

    redshift = 1.0/a_exp - 1


    return redshift, outputlist







#==========================================
def find_most_massives(srcdir, nhalos):
#==========================================

    outputfiles = listdir(srcdir)
    halofiles = []
    for h in outputfiles:
        if 'halo_' in h:
            halofiles.append(h)

    if len(halofiles) == 0:
        print("Found no halo_XXXXX.txtYYYYY files in ", srcdir)


    halos = np.array([], dtype='int')
    masses = np.array([], dtype='float')

    for hfile in halofiles:
        newhalo, newmass = np.loadtxt(srcdir+'/'+hfile, usecols=([0, 6]), dtype='float', skiprows=1, unpack=True)
        halos = np.concatenate((halos, np.atleast_1d(newhalo.astype('int'))))
        masses = np.concatenate((masses, np.atleast_1d(newmass)))


    sortind = masses.argsort()

    # print mass as well
    #  print('{0:15}{1:15}'.format('halo', 'mass'))
    #  print('{0:15}{1:15}'.format('--------------', '--------------'))
    #  for i in range(nhalos):
    #      print('{0:<15d}{1:<15.6E}'.format(halos[sortind[-i-1]], masses[sortind[-i-1]]))

    # print just clump ID
    for i in range(nhalos):
        print('{0:<15d}'.format(halos[sortind[-i-1]]))







#==================================
if __name__ == '__main__':
#==================================

    nhalos, srcdir = parse_args()
    if srcdir == "None":
        # find z = 0 directory
        z, outputdirs = read_redshifts()
        z0ind = np.argmin(abs(z))
        srcdir = outputdirs[z0ind]
    else:
        if not path.exists(srcdir):
            print("Can't find directory", srcdir)
            quit()
    

    find_most_massives(srcdir, nhalos)
