#!/usr/bin/python2

#===============================================
# A simple quick plot of all particles
# in a output_XXXXX directory.
# Usage: quickpartplot.py <output_nr>
# output_nr: number of output, not zero padded,
# not full name
#===============================================





#========================================
def read_particle_data(srcdir, ncpu):
#========================================
    """
    Reads in the particle data from directory srcdir.
    NOTE: requires part_XXXXX.outYYYYY and unbinding_XXXXX.outYYYYY files

    Parameters:
        srcdir:     String of directory where to read data from
        ncpu:       ncpu parameter

    returns:
        x,y:      numpy arrays of particle positions
    """

    import numpy as np
    from os import listdir
    import fortranfile as ff

    srcdirlist = listdir(srcdir)

    if 'part_'+srcdir[-5:]+'.out00001' not in srcdirlist:
        print "Couldn't find unbinding_"+srcdir[-5:]+".out00001 in", srcdir
        print "To plot particles, I require the particle output."
        quit()





    #-----------------------
    # First read headers
    #-----------------------
    nparts = np.zeros(ncpu, dtype='int')
    partfiles = [0]*ncpu

    for cpu in range(ncpu):
        srcfile = srcdir+'/part_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)
        partfiles[cpu] = ff.FortranFile(srcfile)

        ncpu_junk = partfiles[cpu].readInts()
        ndim = partfiles[cpu].readInts()
        nparts[cpu] = partfiles[cpu].readInts()
        localseed = partfiles[cpu].readInts()
        nstar_tot = partfiles[cpu].readInts()
        mstar_tot = partfiles[cpu].readReals('d')
        mstar_lost = partfiles[cpu].readReals('d')
        nsink = partfiles[cpu].readInts()

        del ncpu_junk, ndim, localseed, nstar_tot, mstar_tot, mstar_lost, nsink



    #-------------------
    # Allocate arrays
    #-------------------
    nparttot = nparts.sum()
    x = np.zeros(nparttot, dtype='float')
    y = np.zeros(nparttot, dtype='float')


    #----------------------
    # Read particle data
    #----------------------

    start_ind = np.zeros(ncpu, dtype='int')
    for cpu in range(ncpu-1):
        start_ind[cpu+1] = nparts[cpu] + start_ind[cpu]

    for cpu in range(ncpu):
        x[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')
        y[start_ind[cpu]:start_ind[cpu]+nparts[cpu]] = partfiles[cpu].readReals('d')

        srcfile = srcdir+'/part_'+srcdir[-5:]+'.out'+str(cpu+1).zfill(5)


    return x, y



#=======================
def read_cmdlineargs():
#=======================
    """
    Reads cmdline args correctly.
    """
    from sys import argv
    from os.path import exists

    def abort(msg=''):
        print "I expect exactly 1 argument: Number of output directory for which to plot particles for"
        print msg
        quit()

    if len(argv) != 2:
        abort()

    try:
        dirnr = int(argv[1])
    except ValueError :
        abort('Argument was not an integer')


    srcdir = 'output_'+str(dirnr).zfill(5)
    if not exists(srcdir):
        print "Directory", srcdir, "not found."
        quit()

    return srcdir


#============================
def read_infofile(srcdir):
#============================
    """
    Reads ncpu parameter.
    """
    infofile = srcdir+'/'+'info_'+srcdir[-5:]+'.txt'
    f = open(infofile, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()
    f.close()

    ncpu = int(line[-1])

    return ncpu






#=============================
def plot_particles(x,y,srcdir):
#=============================
    """
    Plots the particles.
    Requires fast-histogram package:
    https://github.com/astrofrog/fast-histogram

    parameters:
        x,y :   numpy arrays of x and y particle positions

    returns:
        nothing
    """

    from fast_histogram import histogram2d
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

    print "Creating figure"

    # imshow wants y-axis as first index
    xmax = x.max()
    xmin = x.min()
    ymax = y.max()
    ymin = y.min()
    if (xmax < 1) and (ymax < 1):
        r = [[0,1],[0,1]]
    else:
        ma = max(xmax, ymax)
        ma = ((int(ma)+1)/10)*10.0
        mi = min(xmin, ymin)
        mi = (int(mi)/10)*10.0
        r = [[mi,ma],[mi,ma]]

    if x.shape[0]>1000000:
        bins = [1000,1000]
    else:
        bins = [200,200]

    hist=histogram2d(y,x, range=r, bins=bins)
    minval = hist[hist>0.0].min()
    hist[hist==0] = minval*5e-1


    fig = plt.figure(figsize=(10,10), dpi=300)
    ax=fig.add_subplot(111)

    im = ax.imshow(hist,
        interpolation='kaiser',
        cmap='inferno',
        origin='lower',
        extent=(r[0][0],r[0][1],r[1][0], r[1][1]),
        norm=LogNorm()
    )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.15)
    fig.colorbar(im, cax=cax)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title("Particle Projection "+srcdir)


    fig.tight_layout()


    print "saving figure"
    plt.savefig('particleplot_'+srcdir[-5:]+'.png')


    return



#===============
def main():
#===============

    srcdir = read_cmdlineargs()
    ncpu = read_infofile(srcdir)

    x, y = read_particle_data(srcdir, ncpu)

    plot_particles(x,y,srcdir)

    #  import subprocess
    #  cmd='eog particleplot_'+srcdir[-5:]+'.png'
    #  exitcode=subprocess.call(cmd, shell=True)




#============================
if __name__=='__main__':
#============================
    main()
