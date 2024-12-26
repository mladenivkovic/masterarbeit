#!/usr/bin/python2

#===========================================
# Plots galaxies from output made by
# get_smf.f03, which is compiled and
# called by get_correlation.sh
# Make sure subroutine write_galaxies is
# called within get_smf.f03
#===========================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import fortranfile as ff


#==========================
def get_outdirlist():
#==========================
    from os import getcwd
    from os import listdir

    workdir = getcwd()
    filelist = listdir(workdir)

    outputlist = []
    outdirtemp = 'output_'
    for filename in filelist:
        if filename.startswith(outdirtemp):
            # first check whether directory contains galaxy files
            txtfilelist = listdir(filename)
            if ('density_image_galaxies-all.dat' in txtfilelist) and ('density_image_galaxies-sub.dat' in txtfilelist):
                outputlist.append(filename)

    if len(outputlist)<1:
        print "I didn't find any output_XXXXX directories in current working directory."
        print "Or none that contain an density_image-all.txt or density_image-sub.txt file."
        print "Are you in the correct workdir?"
        quit()

    outputlist.sort()
    noutput = len(outputlist)
    print "First output containing galaxies: ", outputlist[0]

    a_exp = np.zeros(noutput, dtype='float')

    #-----------------------------
    # Read info files
    #-----------------------------

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

            infofile.close()

        except IOError: # If file doesn't exist
            print "Didn't find any info data in ", srcdir
            break

    redshift = np.round(np.abs(1.0/a_exp - 1), 3)

    return outputlist, redshift



#==========================
def plot_dir(srcdir, z):
#==========================
    """
    Plots image of directory srcdir
    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    nr = srcdir[-5:]

    print "Plotting ", srcdir

    #  for i, j in enumerate([('all', 'including orphans'), ('sub', 'excluding orphans')]):
    #
    #      case = j[0]
    #      title = j[1]
    #
    #      fname = srcdir+'/density_image-'+case+'.txt'
    #      image = np.loadtxt(fname, dtype='float')
    #
    #      print "min:\t", np.log10(image[image>0].min()), "\t max:\t", np.log10(image.max())
    #

    fig = plt.figure(figsize=(18,12), dpi=300)
    ax1 = fig.add_subplot(2,3,1)
    ax2 = fig.add_subplot(2,3,2)
    ax3 = fig.add_subplot(2,3,3)
    ax4 = fig.add_subplot(2,3,4)
    ax5 = fig.add_subplot(2,3,5)
    ax6 = fig.add_subplot(2,3,6)
    upper_row = [ax1, ax2, ax3]
    lower_row = [ax4, ax5, ax6]
    rows = [upper_row, lower_row]

    imagedata = []

    # first read in everything so we have same min/max values for each plot
    minval = 1e30
    maxval = -1.

    for i, j in enumerate([('all', 'including orphans'), ('sub', 'excluding orphans')]):

        case = j[0]

        fname = srcdir+'/density_image_galaxies-'+case+'.dat'
        f = ff.FortranFile(fname)
        nc = np.asscalar(f.readInts())
        aexp = np.asscalar(f.readReals('d'))
        unit_l = np.asscalar(f.readReals('d'))
        xmax = np.asscalar(f.readReals('d'))
        xmin = np.asscalar(f.readReals('d'))
        ymax = np.asscalar(f.readReals('d'))
        ymin = np.asscalar(f.readReals('d'))
        zmax = np.asscalar(f.readReals('d'))
        zmin = np.asscalar(f.readReals('d'))
        imageXY = f.readReals('d')
        imageYZ = f.readReals('d')
        imageXZ = f.readReals('d')

        imageXY = imageXY.reshape((nc,nc))
        imageYZ = imageYZ.reshape((nc,nc))
        imageXZ = imageXZ.reshape((nc,nc))

        imagedata.append([imageXY, imageYZ, imageXZ])

        for im in [imageXY, imageYZ, imageXZ]:
            minval = min(im[im>0].min(), minval)
            maxval = max(im.max(), maxval)

    minval *= 1e-3

    for i, j in enumerate([('all', 'including orphans'), ('sub', 'excluding orphans')]):

        case = j[0]
        title = j[1]

        imageXY, imageYZ, imageXZ = imagedata[i]

        imageXY[imageXY==0] = minval
        imageYZ[imageYZ==0] = minval
        imageXZ[imageXZ==0] = minval


        # Plot image
        for ax in rows[i]:
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        ax1, ax2, ax3 = rows[i]
        ax1.imshow(imageXY,
            origin='lower',
            extent=(0,1,0,1),
            norm=colors.LogNorm(),
            interpolation='gaussian',
            cmap='magma',
            vmin=minval, vmax=maxval
            )
        ax1.set_xlabel("x")
        ax1.set_ylabel("y")

        ax2.imshow(imageYZ,
            origin='lower',
            extent=(0,1,0,1),
            norm=colors.LogNorm(),
            interpolation='gaussian',
            cmap='magma',
            vmin=minval, vmax=maxval
            )
        ax2.set_xlabel("y")
        ax2.set_ylabel("z")
        ax2.set_title(title)

        ax3.imshow(imageXZ, 
            origin='lower',
            extent=(0,1,0,1),
            norm=colors.LogNorm(),
            interpolation='gaussian',
            cmap='magma',
            vmin=minval, vmax=maxval
            )
        ax3.set_xlabel("x")
        ax3.set_ylabel("z")

    plt.subplots_adjust( 
                        left = 0.025,   # the left side of the subplots of the figure
                        right = 0.98,    # the right side of the subplots of the figure
                        bottom = 0.05,   # the bottom of the subplots of the figure
                        top = 0.95,      # the top of the subplots of the figure
                        wspace = 0.1,   # the amount of width reserved for space between subplots,
                                      # expressed as a fraction of the average axis width
                        hspace = 0.1,   # the amount of height reserved for space between subplots,
                                      # expressed as a fraction of the average axis height
                        )
    plt.figtext(.02, .03, 'z={0:5.3f}'.format(z))

    outfname = 'galaxy_density_projection-'+nr+'.png'
    plt.savefig(outfname, dpi=300)
    plt.close()




#===================
def main():
#===================

    dirlist, redshift = get_outdirlist()

    for i in range(len(dirlist)):
        plot_dir(dirlist[i], redshift[i])





#==============================
if __name__ == "__main__":
#==============================
    main()

