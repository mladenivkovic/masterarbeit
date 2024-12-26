#!/usr/bin/python2

#------------------------------------------
# Module containing variables and functions
# used for various evaluations.
#------------------------------------------

from __future__ import print_function
import numpy as np
import warnings


#==============================================================================
#============== GLOBAL PARAMS
#==============================================================================


#===========================
class global_vars:
#===========================
    """
    A class that contains globally used variables.
    """

    noutput = 0
    lastdir = 'output_00001'
    ncpu = 1
    workdir = '~/'
    outdirtemp = 'output_'
    dftemp = 'galaxies_'
    ozi = 0


    outputnrs = []
    a_exp = []
    redshift = []
    unit_l = []
    unit_m = []
    unit_t = []
    unit_dens = []
    h = 0


    galaxy_clumps = 0
    galaxy_masses = 0
    galaxy_pos = 0
    halos = 0
    halo_masses = 0





#=============================
class physics_constants:
#=============================
    """
    A class to store physical constants
    """

    Mpc = 3.086e24 # cm



g = global_vars()
ph = physics_constants()





#==============================================================================
#============ READING FILES
#==============================================================================


#==========================
def get_output_info():
#==========================
    """
    Read in the output info based on the files in the current
    working directory.
    Reads in last directory, ncpu, noutputs, workdir.

    parameters:     none
    returns:        nothing.
    """

    from os import getcwd
    from os import listdir

    global g

    g.workdir = getcwd()
    filelist = listdir(g.workdir)

    outputlist = []
    for filename in filelist:
        if g.outdirtemp in filename:
            outputlist.append(filename)


    if len(outputlist)<1:
        print("I didn't find any output_XXXXX directories in current working directory.")
        print("Are you in the correct workdir?")
        quit()

    outputlist.sort()
    g.lastdir = outputlist[-1]
    g.noutput = len(outputlist)

    # read ncpu from infofile in last output directory
    infofile = g.lastdir+'/'+'info_'+g.lastdir[-5:]+'.txt'
    f = open(infofile, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()

    g.ncpu = int(line[-1])
    f.close()

    return






#========================================
def read_info_files():
#========================================
    """
    Reads info files.
    """

    import numpy as np


    global g

    g.outputnrs = np.array(range(g.noutput, 0, -1))
    g.a_exp = np.zeros(g.noutput)
    g.unit_l = np.zeros(g.noutput)
    g.unit_m = np.zeros(g.noutput)
    g.unit_t = np.zeros(g.noutput)
    g.unit_dens = np.zeros(g.noutput)

    fn_temp = 'info_'

    for out in range(g.noutput):
        dirnr = str(g.noutput-out).zfill(5)
        dirnm = g.outdirtemp+dirnr+'/'
        fn = dirnm+fn_temp+dirnr+".txt"

        try:
            infofile = open(fn)
            for i in range(9):
                infofile.readline() # skip first 9 lines

            # get expansion factor
            aline = infofile.readline()
            astring, equal, aval = aline.partition("=")
            afloat = float(aval)
            g.a_exp[out] = afloat

            # get H0 factor
            hline = infofile.readline()
            hstring, equal, hval = hline.partition("=")
            hfloat = float(hval)
            g.h = hfloat/100.0

            for i in range(4):
                infofile.readline() # skip 5 lines

            # get unit_l
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            g.unit_l[out] = unitfloat

            # get unit_dens
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            g.unit_dens[out] = unitfloat

            # get unit_t
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            g.unit_t[out] = unitfloat

        except IOError: # If file doesn't exist
            print("Didn't find any info data in ", srcdir)
            break


    g.redshift = 1.0/g.a_exp - 1
    g.ozi = np.argmin(np.absolute(g.redshift))
    print("z=0 snapshot: ", g.outputnrs[g.ozi])

    return






#================================================
def read_data(read_pos=True, logmasses=False):
#================================================
    """
    Reads galaxy and halo data
    Parameters:
        read_pos:   Boolean whether to read galaxy positions
    Todo: dox
    """
    global g

    outputnr = g.outputnrs[g.ozi]
    on = str(outputnr).zfill(5)
    srcdir = g.outdirtemp+on+'/'
    fn = g.dftemp+on+'.txt' # gives galaxies_XXXXX.txt


    #--------------------------
    # Read galaxy data
    #--------------------------

    try:
        fg = [srcdir+fn+str(cpu+1).zfill(5) for cpu in range(g.ncpu)]

        clumps = [np.zeros((1)) for i in range(g.ncpu)]
        masses = [np.zeros((1)) for i in range(g.ncpu)]
        if read_pos :
            pos = [np.zeros((1,3)) for i in range(g.ncpu)]


        i = 0
        for f in fg:
            with warnings.catch_warnings():
                warnings.filterwarnings('error') # treat warnings as errors so I can catch them
                try:
                    clumps[i] = np.atleast_1d(np.loadtxt(f, dtype='int', skiprows=1, usecols=([0])))
                    masses[i] = np.atleast_1d(np.loadtxt(f, dtype='float', skiprows=1, usecols=([1])))
                    if (read_pos):
                        pos[i] = np.loadtxt(f, dtype='float', skiprows=1, usecols=([2, 3, 4]), unpack=False)
                    i += 1
                except Warning:
                    continue
        if i > 0:
            g.galaxy_clumps = np.concatenate(clumps[:i], axis=0)
            if logmasses:
                g.galaxy_masses = np.log10(np.concatenate(masses[:i], axis=0))
            else:
                g.galaxy_masses = np.concatenate(masses[:i], axis=0)
            if (read_pos):
                g.galaxy_pos = np.concatenate(pos[:i], axis=0)


    except IOError: # If file doesn't exist
        print("Didn't find any galaxy data in", srcdir)
        quit()



    #--------------------------
    # Read halo data
    #--------------------------

    fn = 'halo_'+on+'.txt' # gives halo_XXXXX.txt
    try:
        fg = [srcdir+fn+str(cpu+1).zfill(5) for cpu in range(g.ncpu)]

        datalist = [np.zeros((1)) for i in range(g.ncpu)]

        i = 0
        for f in fg:
            with warnings.catch_warnings():
                warnings.filterwarnings('error') # treat warnings as errors so I can catch them
                try:
                    datalist[i] = np.atleast_1d(np.loadtxt(f, dtype='int', skiprows=1, usecols=([0])))
                    i += 1
                except Warning:
                    continue
        if i > 0:

            g.halos = np.concatenate(datalist[:i], axis=0)


    except IOError: # If file doesn't exist
        print("Didn't find any halo data in", srcdir)
        quit()



#      #--------------------------
    #  # Read mergertree data
    #  #--------------------------
    #
    #  fn = 'mergertree_'+on+'.txt' # gives halo_XXXXX.txt
    #  try:
    #      fg = [srcdir+fn+str(cpu+1).zfill(5) for cpu in range(ncpu)]
    #
    #      datalist = [np.zeros((1,3)) for i in range(ncpu)]
    #
    #      i = 0
    #      for f in fg:
    #          with warnings.catch_warnings():
    #              warnings.filterwarnings('error') # treat warnings as errors so I can catch them
    #              try:
    #                  datalist[i] = np.atleast_1d(np.loadtxt(f, dtype='float', skiprows=1, usecols=([0, 2])))
    #                  i += 1
    #              except Warning:
    #                  continue
    #      if i > 0:
    #          fulldata = np.concatenate(datalist[:i], axis=0)
    #          mtree_halos = fulldata[0][:]
    #          masses = fulldata[0][:]
    #
    #
    #  except IOError: # If file doesn't exist
    #      print("Didn't find any halo data in", srcdir)
    #      quit()
#



    #  #--------------------------
    #  # transform data
    #  #--------------------------
    #  Mpc = 3.086e24 # cm
    #  M_Sol = 1.98855e33 # g
    #
    #  vol = (unit_l[oind]/Mpc)**3
    #  unit_m = unit_dens[oind] * unit_l[oind]**3/M_Sol
    #
    #  #  masses *= unit_m
    #  masses = np.log10(masses)


    return





#========================================
def read_particle_data():
#========================================
    """
    Reads in the particle data from directory srcdir.
    NOTE: requires part_XXXXX.outYYYYY and unbinding_XXXXX.outYYYYY files

    Parameters:
        None

    returns:
        x,y,z:      numpy arrays of particle positions
        clumpid:    numpy arrays of particle clump IDs
    """

    from os import listdir
    import fortranfile as ff

    srcdir = g.outdirtemp+str(g.outputnrs[g.ozi]).zfill(5)
    srcdirlist = listdir(srcdir)

    if 'unbinding_'+srcdir[-5:]+'.out00001' not in srcdirlist:
        print("Couldn't find unbinding_"+srcdir[-5:]+".out00001 in", srcdir)
        print("To plot particles, I require the unbinding output.")
        quit()






    #-----------------------
    # First read headers
    #-----------------------
    nparts = np.zeros(g.ncpu, dtype='int')
    partfiles = [0]*g.ncpu

    for cpu in range(g.ncpu):
        srcfile = srcdir+'/part_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)
        partfiles[cpu] = ff.FortranFile(srcfile)

        ncpu = partfiles[cpu].readInts()
        ndim = partfiles[cpu].readInts()
        nparts[cpu] = partfiles[cpu].readInts()
        localseed = partfiles[cpu].readInts()
        nstar_tot = partfiles[cpu].readInts()
        mstar_tot = partfiles[cpu].readReals('d')
        mstar_lost = partfiles[cpu].readReals('d')
        nsink = partfiles[cpu].readInts()

        del ncpu, ndim, localseed, nstar_tot, mstar_tot, mstar_lost, nsink



    #-------------------
    # Allocate arrays
    #-------------------
    nparttot = nparts.sum()
    x = np.zeros(nparttot, dtype='float')
    y = np.zeros(nparttot, dtype='float')
    z = np.zeros(nparttot, dtype='float')
    clumpid = np.zeros(nparttot, dtype='int')


    #----------------------
    # Read particle data
    #----------------------

    start_ind = np.zeros(g.ncpu, dtype='int')
    for cpu in range(g.ncpu-1):
        start_ind[cpu+1] = nparts[cpu] + start_ind[cpu]

    for cpu in range(g.ncpu):
        x[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')
        y[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')
        z[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')

        srcfile = srcdir+'/part_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)
        unbfile = srcdir+'/unbinding_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)
        unbffile = ff.FortranFile(unbfile)

        clumpid[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = unbffile.readInts()


    clumpid = np.absolute(clumpid)


    return x, y, z, clumpid










#==============================================================================
#=============== CALCULATING
#==============================================================================


#=======================================
def get_density_field(nc, which='all'):
#=======================================
    """
    Calculates the density field.
    nc : number of cells for field per dimension
    which:  all:    take mains, subhaloes and orphans
            main:   main only
            sub:    main and sub, no orphans
    """

    global g

    density_field = np.zeros((nc, nc, nc), dtype='float')

    cellvolume = (g.unit_l[g.ozi]/ph.Mpc/nc)**3
    if which=='all':
        density_field, bins = np.histogramdd(g.galaxy_pos, bins=(nc,nc,nc), range=((0,1), (0,1), (0,1)), weights=g.galaxy_masses)

    elif which=='main':
        i = 0
        pos = np.zeros(g.galaxy_pos.shape, dtype='float')
        mass  = np.zeros(g.galaxy_masses.shape, dtype='float')
        for j, cid in enumerate(g.galaxy_clumps):
            if cid in g.halos:
                pos[i] = g.galaxy_pos[j]
                mass[i] = g.galaxy_masses[j]
                i += 1

        density_field, bins = np.histogramdd(pos[:i+1], bins=(nc,nc,nc), range=((0,1), (0,1), (0,1)), weights=mass[:i+1])

    elif which == 'sub':
        i = 0
        pos = np.zeros(g.galaxy_pos.shape, dtype='float')
        mass  = np.zeros(g.galaxy_masses.shape, dtype='float')
        for j, cid in enumerate(g.galaxy_clumps):
            if cid > 0:
                pos[i] = g.galaxy_pos[j]
                mass[i] = g.galaxy_masses[j]
                i += 1

        density_field, bins = np.histogramdd(pos[:i+1], bins=(nc,nc,nc), range=((0,1), (0,1), (0,1)), weights=mass[:i+1])

    density_field /= cellvolume




    return density_field



