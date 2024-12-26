#!/usr/bin/env python3

#======================================================
# Compute and plot satellite profile for a halo
# usage:
#   satellite_profile.py output_XXXXX halo-id
#======================================================


import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

errormsg='''
Compute and plot satellite profile for a halo.
Usage:
    satellite_profile.py output_XXXXX halo-id
or:
    satellite_profile.py --fb path/to/file.pkl

The program, when executed once, will dump a .pkl
file with all the necessary data to make a plot
in the current workdir.
'''


# Plot parameters
params = {
    'axes.labelsize': 12,
    'axes.titlesize': 12*1.1,
    'font.size': 12,
    'font.family': 'serif',
    'font.serif': 'Computer Modern',
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.subplot.left'    : 0.05,
    'figure.subplot.right'   : 0.97,
    'figure.subplot.bottom'  : 0.10,
    'figure.subplot.top'     : 0.92,
    'figure.subplot.wspace'  : 0.15,
    'figure.subplot.hspace'  : 0.15,
    'figure.dpi' : 300,
    'lines.markersize' : 6,
    'lines.linewidth' : 2.
    }


mpl.rcParams.update(params)




#======================
class global_params:
#======================
    """
    An object to store all global parameters in, so you can only
    pass 1 object to functions that need it.
    """


    #======================
    def __init__(self):
    #======================
        """
        Initialises object.
        """

        self.halo           = 0     # the root halo for which to plot for
        self.dirname        = ''    # output directory to work for
        self.outputnr       = 0     # XXXXX in output_XXXXX
        self.ncpu           = 0     # how many processors were used for simulation run
        self.outputfilename = ''    # filename for finished tree image
        self.prefix         = ''    # filename prefix
        self.workdir        = ''    # current working directory
        self.from_backup    = False # whether we're using a backup pickle file for the data
        self.picklefile     = ''    # path to pickle file


        # Define a long list of colors for the branches
        self.colorlist=[
                'black',
                'red',
                'green',
                'gold',
                'cornflowerblue',
                'lime',
                'magenta',
                'orange',
                'mediumpurple',
                'deeppink',
                'lightgreen',
                'saddlebrown',
                'orchid',
                'mediumseagreen']

        return





    #=============================
    def read_cmdlineargs(self):
    #=============================
        """
        Reads in the command line arguments and stores them in the
        global_params object.
        """
        from sys import argv
        import os.path
        from os import getcwd
        from os import listdir

        nargs = len(argv)

        if nargs != 3:
            print("Need exactly 2 args: outputdir and halo ID, or -fb and .pkl file")
            quit()

        if '-h' in argv or '--help' in argv:
            print(errormsg)
            quit()

        if '-fb' in argv or '--fb' in argv:
            self.from_backup = True
            if argv[1] == '-fb' or argv[1] == '--fb':
                self.picklefile = argv[2]
            else:
                self.picklefile = argv[1]

            if not os.path.isfile(self.picklefile):
                print("file:", self.picklefile, "could not be found.")
                quit()
            return

        self.dirname = argv[1]
        self.dirname = self.dirname.rstrip('/')

        if not os.path.isdir(self.dirname):
            print(self.dirname, "is not a directory. Try again.")
            quit()

        try:
            self.halo = int(argv[2])
        except ValueError:
            print("I didn't recognize the argument '", arg, "'")
            print("use satellite_profile.py -h or --help to print help message.")
            quit()

        if self.halo <= 0:
            print("No or wrong halo given. Halo ID must be > 0")
            quit()



        # read in other metadata
        self.workdir = getcwd()
        self.outputnr = self.dirname[-5:]
        try:
            nr = int(self.outputnr)
        except ValueError:
            print("Something went wrong when fetching the output number from", self.dirname[-5:])
            quit()


        # read ncpu from infofile
        infofile = self.dirname+'/'+'info_'+self.dirname[-5:]+'.txt'
        f = open(infofile, 'r')
        ncpuline = f.readline()
        line = ncpuline.split()

        self.ncpu = int(line[-1])
        f.close()


        print("Working for dir", self.dirname, "and halo", self.halo, 'with ncpu =', self.ncpu)


        return




#======================================
def get_clump_data(params):
#======================================
    """
    reads in clumpfinder data from files.
    params: global_params object

    returns:
        data_int: numpy array of integer data
            data_int[0]: clump ID
            data_int[1]: clump level
            data_int[2]: ID of parent clump
    """


    print("Reading in clump data.")


    datalist = [np.zeros(3) for cpu in range(params.ncpu)]

    # Loop over files
    for cpu in range(params.ncpu):
        inputfile = params.dirname+'/clump_'+params.outputnr.zfill(5)+'.txt'+str(cpu+1).zfill(5)

        # get clump ids, parents and levels
        datalist[cpu] = np.loadtxt(inputfile, dtype='int', skiprows=1, usecols=[0,1,2])


    data_int = np.concatenate(datalist, axis=0)

    return data_int







#==============================================
def get_all_children(params, clumpdata):
#==============================================
    """
    Find children recursively for given halo.

    params: global_params object
    clumpdata: result of get_clump_data

    returns:
        children:       sorted list of all children of halo
    """

    # identify all children and their levels recursively for given halo
    children = []
    i = 0
    while(i < clumpdata.shape[0]):
        try:
            if (clumpdata[i, 2] in children or clumpdata[i,2]==params.halo): #if parent is halo or one of the children
                if (clumpdata[i,0] != params.halo): #if clump itself isn't halo
                    if (clumpdata[i,0] not in children):
                        children.append(clumpdata[i,0])
                        i = -1                # start again from the beginning to find hierarchical buildup

        except IndexError: # in case there is only one clump in simulation
            children = []
            break

        i+=1

    return sorted(children)








#========================================
def read_particle_data(params, cpu):
#========================================
    """
    Reads in the particle data from file written by cpu cpu+1.
    NOTE: requires part_XXXXX.outYYYYY and unbinding_XXXXX.outYYYYY files

    Parameters:
        params:     global_params() object
        cpu:        for which cpu to work

    returns:
        partID:     particle ID
        clumpid:    numpy arrays of particle clump IDs
    """

    from scipy import io as ff



    #-----------------------
    # First read headers
    #-----------------------
    srcfile = params.dirname+'/part_'+str(params.outputnr).zfill(5)+'.out'+str(cpu+1).zfill(5)
    pyfile =  ff.FortranFile(srcfile)

    ncpu = pyfile.read_ints()
    ndim = pyfile.read_ints()
    nparts = pyfile.read_ints()
    localseed = pyfile.read_ints()
    nstar_tot = pyfile.read_ints()
    mstar_tot = pyfile.read_reals('d')
    mstar_lost = pyfile.read_reals('d')
    nsink = pyfile.read_ints()




    #----------------------
    # Read particle data
    #----------------------

    x = pyfile.read_reals('d')
    y = pyfile.read_reals('d')
    z = pyfile.read_reals('d')
    vx = pyfile.read_reals('d')
    vy = pyfile.read_reals('d')
    vz = pyfile.read_reals('d')
    mass = pyfile.read_reals('d')
    partID = pyfile.read_ints()

    del x, y, z, vx, vy, vz, mass

    unbfile = params.dirname+'/unbinding_'+str(params.outputnr).zfill(5)+'.out'+str(cpu+1).zfill(5)
    unbffile = ff.FortranFile(unbfile)

    clumpid = unbffile.read_ints()
    clumpid = np.absolute(clumpid)

    return partID, clumpid









#======================================
def read_galaxy_data(params):
#======================================
    """
    reads in galaxy data as written by the mergertree patch.
    NOTE: requires galaxies_XXXXX.txtYYYYY files.

    parameters:
        params:         global_params() object

    returns:
        galdata = [gal_halo, gal_mass, gal_pos]
        orphdata = [orph_partid, orph_mass, orph_pos]
            lists of np.arrays, galaxy clump, stellar mass, and positions,
            and orphan associated particle, stellar mass, position
        centerpos: np.array(3), position of central galaxy
    """

    import warnings
    import gc
    from os import listdir

    print("Reading in galaxy data.")

    # define new datatype for galaxy output
    galtype = np.dtype([  ('clump', 'i4'),
                        ('mass', 'f8'),
                        ('x', 'f8'),
                        ('y', 'f8'),
                        ('z', 'f8'),
                        ('partID', 'i4')
                        ])



    datalist = [None for i in range(params.ncpu)]

    for cpu in range(params.ncpu):
        srcfile = params.dirname+'/galaxies_'+str(params.outputnr).zfill(5)+'.txt'+str(cpu+1).zfill(5)

        datalist[cpu] = np.atleast_1d(np.genfromtxt(srcfile, dtype=galtype, skip_header=1))


    # count the number of galaxies in each file
    ngals = [0 for cpu in range(params.ncpu)]

    for cpu in range(params.ncpu):
        i = 0
        while datalist[cpu][i]['clump'] > 0:
            i+=1
        ngals[cpu] = i

    ngals_tot = sum(ngals)
    norph_tot = sum( [datalist[cpu].shape[0]-ngals[cpu] for cpu in range(params.ncpu)] )



    gal_pos =       np.zeros((ngals_tot, 3))
    gal_mass =      np.zeros((ngals_tot))
    gal_halo =      np.zeros((ngals_tot), dtype=np.int)
    orph_pos =      np.zeros((norph_tot, 3))
    orph_mass =     np.zeros((norph_tot))
    orph_partid =   np.zeros((norph_tot), dtype=np.int)

    ind_gal = 0
    ind_orph = 0

    for cpu in range(params.ncpu):

        gal_halo[ind_gal:ind_gal+ngals[cpu]] = datalist[cpu][:ngals[cpu]]['clump']
        gal_mass[ind_gal:ind_gal+ngals[cpu]] = datalist[cpu][:ngals[cpu]]['mass']
        gal_pos[ind_gal:ind_gal+ngals[cpu], 0] = datalist[cpu][:ngals[cpu]]['x']
        gal_pos[ind_gal:ind_gal+ngals[cpu], 1] = datalist[cpu][:ngals[cpu]]['y']
        gal_pos[ind_gal:ind_gal+ngals[cpu], 2] = datalist[cpu][:ngals[cpu]]['z']
        ind_gal += ngals[cpu]

        l = datalist[cpu].shape[0]
        orph_mass[ind_orph:ind_orph+l-ngals[cpu]] = datalist[cpu][ngals[cpu]:]['mass']
        orph_pos[ind_orph:ind_orph+l-ngals[cpu], 0] = datalist[cpu][ngals[cpu]:]['x']
        orph_pos[ind_orph:ind_orph+l-ngals[cpu], 1] = datalist[cpu][ngals[cpu]:]['y']
        orph_pos[ind_orph:ind_orph+l-ngals[cpu], 2] = datalist[cpu][ngals[cpu]:]['z']
        orph_partid[ind_orph:ind_orph+l-ngals[cpu]] = datalist[cpu][ngals[cpu]:]['partID']
        ind_orph += l - ngals[cpu]


    galdata = [gal_halo, gal_mass, gal_pos]
    orphdata = [orph_partid, orph_mass, orph_pos]


    # find central galaxy position

    for i, g in enumerate(gal_halo):
        if g == params.halo:
            centerpos = gal_pos[i]
            break


    return galdata, orphdata, centerpos






#==============================================================
def get_subhalo_galaxies(params, children, galaxies):
#==============================================================
    """
    Get a list of masses and positions of satellite galaxies
    for given halo
    """

    masses = np.empty(10000)
    positions = np.empty((10000, 3))

    gid = galaxies[0]
    gm = galaxies[1]
    gpos = galaxies[2]

    ngals = 0
    for i,g in enumerate(gid):
        if g in children:
            masses[ngals] = gm[i]
            positions[ngals] = gpos[i]
            ngals += 1

            if ngals == masses.shape[0]:
                masses.resize(masses.shape[0] + 10000)
                positions.resize((positions.shape[0] + 10000, 3))


    return masses[:ngals], positions[:ngals]







#============================================================
def get_orphan_galaxies(params, children, orphans):
#============================================================
    """
    Get the orphan galaxies that are in this clump.
    """


    print("Sorting out orphans")

    masses = np.empty(10000)
    positions = np.empty((10000, 3))

    oid = orphans[0]
    om = orphans[1]
    opos = orphans[2]

    indsort = np.argsort(oid)
    oid = oid[indsort]
    om = om[indsort]
    opos = opos[indsort]



    ngals = 0

    for cpu in range(params.ncpu):
        # first read in paricle clump ID's and particle IDs
        pids, clumpids = read_particle_data(params, cpu)

        indsort = np.argsort(pids)
        pids = pids[indsort]
        clumpids = clumpids[indsort]

        # now loop over all orphan particles

        i = 0
        j = 0

        while i < oid.shape[0] and j < pids.shape[0]:

            if oid[i] == pids[j]:

                if clumpids[j] > 0:
                    if clumpids[j] in children:
                        masses[ngals] = om[i]
                        positions[ngals] = opos[i]
                        ngals += 1

                        if ngals == masses.shape[0]:
                            masses.resize(masses.shape[0] + 10000)
                            positions.resize((positions.shape[0] + 10000, 3))

                i += 1
                j += 1


            elif oid[i] > pids[j]:
                j += 1
            elif oid[i] < pids[j]:
                i += 1

    return masses[:ngals], positions[:ngals]







#=============================================================================
def dump_backup(params, submass, subpos, orphmass, orphpos, centerpos):
#=============================================================================
    """
    Dump a pickle backup so you can skip the reading and sorting
    if you're only tweaking the plotting
    """

    import pickle
    fname='satellite-profile-output_'+str(params.outputnr).zfill(5)+'-halo-'+str(params.halo)+'.pkl'

    print("Dumping pickle")

    f = open(fname, 'wb')
    pickle.dump([params, submass, subpos, orphmass, orphpos, centerpos], f)
    f.close()

    return



#=============================
def read_backup(params):
#=============================
    """
    Read in the pickle backup
    returns: params, submass, subpos, orphmass, orphpos, centerpos
    """

    import pickle
    f = open(params.picklefile, 'rb')
    return pickle.load(f)












#==================================================================
def plot_profiles(params, gmass, gpos, omass, opos, centerpos):
#==================================================================
    """
    Plot number density of satellites without orphans
    gmass, gpos: np.arrays of galaxy masses and positions
    omass, opos: np.arrays of orphan masses and positions
    centerpos: positions of central galaxy
    """

    nbins = 100


    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    ax1.set_title("Number of Satellites")
    ax2.set_title("Cumulative Number of Satellites")


    print(gpos)
    print(opos)
    allpos = np.concatenate((gpos, opos))
    allmass = np.concatenate((gmass, omass), axis=0)

    rgal = np.empty(gmass.shape)
    for i, pos in enumerate(gpos):
        x = pos - centerpos
        rgal[i] = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    rall = np.empty(allmass.shape)
    for i, pos in enumerate(allpos):
        x = pos - centerpos
        rall[i] = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)

    rmax = max(rgal.max(), rall.max())
    rgal /= rmax
    rall /= rmax


    mthresh = [0, 1e8, 1e9]
    mthreshlabels = ['No mass threshold', r'$M > 10^8 M_{\odot}$', r'$M > 10^9 M_{\odot}$']

    for i, t in enumerate(mthresh):

        hist_gal, bin_edges = np.histogram(rgal[gmass>t], bins = nbins, range=(0, 1))
        hist_all, bin_edges = np.histogram(rall[allmass>t], bins = nbins, range=(0, 1))
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

        ax1.semilogy(bin_centers, hist_gal,
            label=mthreshlabels[i], c=params.colorlist[i], ls='-')
        ax1.semilogy(bin_centers, hist_all,
            label=mthreshlabels[i]+', including orphans', c=params.colorlist[i], ls='--')

        # get cumulative histograms
        for j in range(1, hist_gal.shape[0]):
            hist_gal[j] += hist_gal[j-1]
        for j in range(1, hist_all.shape[0]):
            hist_all[j] += hist_all[j-1]

        ax2.semilogy(bin_centers, hist_gal,
            label=mthreshlabels[i], c=params.colorlist[i], ls='-')
        ax2.semilogy(bin_centers, hist_all,
            label=mthreshlabels[i]+', including orphans', c=params.colorlist[i], ls='--')


    ax1.legend()
    ax2.legend()
    ax1.grid()
    ax2.grid()
    ax1.set_xlabel(r'$r/r_{max}$')
    ax2.set_xlabel(r'$r/r_{max}$')


    plt.savefig('satellite-profile-output_'+str(params.outputnr).zfill(5)+'-halo-'+str(params.halo)+'.png', dpi=300)
    print("finished saving fig")

    return






#=======================
def main():
#=======================

    params = global_params()
    params.read_cmdlineargs()


    if not params.from_backup:
        clumpdata = get_clump_data(params)
        children = get_all_children(params, clumpdata)
        galaxies, orphans, centerpos = read_galaxy_data(params)

        submass, subpos = get_subhalo_galaxies(params, children, galaxies)
        orphmass, orphpos = get_orphan_galaxies(params, children, orphans)

        dump_backup(params, submass, subpos, orphmass, orphpos, centerpos)

    else:
        params, submass, subpos, orphmass, orphpos, centerpos = read_backup(params)

    plot_profiles(params, submass, subpos, orphmass, orphpos, centerpos)



if __name__ == '__main__':
    main()
