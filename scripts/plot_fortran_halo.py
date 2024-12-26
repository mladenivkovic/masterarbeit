#!/usr/bin/python2
# need python2 here for fortranfile module

#===========================================
# Plots halos from output made by
# get_halostats.f03, which is compiled and
# called by get_halostats.sh
# Make sure subroutine write_halo_image() is
# called within get_halostats.f03
#===========================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import fortranfile as ff

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif',
    'serif':['Computer Modern Roman'],
    'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)


#  horizontal=False
horizontal=True
plot_galaxies = False
use_internal_units = True
use_partcenter = True
#  use_internal_units = False


#==========================================
def plot_halo_and_galaxies(dirnr, halonr):
#==========================================
    """
    Plots pdf image of halo and galaxies.
    Needs output_XXXXX number and halo number.
    """
    # NOTE: THIS IS FOR OLD VERSION (PRE 04.2021)

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    nr = str(dirnr).zfill(5)

    #  print("Plotting ", halonr, "from output_"+nr)





    fname = "output_"+nr+'/density_image_halo-'+str(halonr)+'.dat'
    f = ff.FortranFile(fname)
    nc = np.asscalar(f.readInts())
    aexp   = np.asscalar(f.readReals('d'))
    unit_l = np.asscalar(f.readReals('d'))
    xmax   = np.asscalar(f.readReals('d'))
    xmin   = np.asscalar(f.readReals('d'))
    ymax   = np.asscalar(f.readReals('d'))
    ymin   = np.asscalar(f.readReals('d'))
    # YOU SHOULD ADAPT THIS TO THE NEW VERSION IF YOU WANT TO REDO THESE PLOTS.
    # SEE THE FUNCTION plot_halo()
    image  = f.readReals('d')
    minval=np.min(image[image>0])
    maxval = np.max(image)
    minval/=100
    image[image==0] = minval
    image = image.reshape((nc,nc))

    # cheat to not have white spaces with lognorm on image
    #  image = np.swapaxes(image, 0, 1)

    # Plot image
    if horizontal:
        fig = plt.figure(figsize=(18,6.2))
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)
    else:
        fig = plt.figure(figsize=(6.2,18))
        ax1 = fig.add_subplot(3,1,1)
        ax2 = fig.add_subplot(3,1,2)
        ax3 = fig.add_subplot(3,1,3)


    #------------------------------
    # Plot Halo Backgrounds
    #------------------------------

    for ax in [ax1, ax2, ax3]:

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.imshow(
            image,
            origin='lower',
            extent=(xmin, xmax, ymin, ymax),
            norm=colors.LogNorm(),
            interpolation='Kaiser',
            #  interpolation='bicubic',
            #  interpolation='gaussian',
            #  cmap='plasma',
            #  cmap='winter',
            #  cmap='bone',
            #  cmap='PiYG',
            #  cmap='PuBu',
            #  cmap='RdBu',
            #  cmap='GnBu',
            cmap='Blues',
            #  cmap='cool',
            #  cmap='inferno',
            #  cmap='magma',
            #  cmap='ocean',
            #  cmap='terrain',
            vmin=minval/100, vmax=maxval
            )

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])

        dx = xmax - xmin
        xs = xmin + 0.05*dx
        xe = xmin + 0.25*dx
        ys = ymin + 0.05*dx
        ye = ymin + 0.07*dx
        xm = 0.5*(xs+xe)

        plt.sca(ax)
        plt.annotate("",
            (xs, ys), (xe, ys),
            arrowprops={
                'arrowstyle':'<->',
                'color':'m',
                'linewidth':'1'
                }
            )
        plt.annotate("{0:3.1f} Mpc".format(dx*unit_l),
            xy=(xm, ye),
            xycoords='data',
            color='m',
            horizontalalignment='center',
            fontsize=15
            )



    #-----------------------
    # Plot Orphan Galaxies
    #-----------------------

    fname = "output_"+nr+'/positions_orphans-'+str(halonr)+'.dat'
    xo, yo = np.loadtxt(fname, dtype='float', skiprows=1, usecols=[0,1], unpack=True)
    ax.scatter(xo, yo,
        facecolor="None",
        edgecolor='orange',
        s=12,
        linewidth=1,
        label='orphan galaxies'
        )




    #-----------------------
    # Plot Clump Galaxies
    #-----------------------

    for ax in [ax2, ax3]:

        fname = "output_"+nr+'/positions_galaxies-'+str(halonr)+'.dat'
        xg, yg = np.loadtxt(fname, dtype='float', skiprows=1, usecols=[0,1], unpack=True)
        ax.scatter(xg[1:], yg[1:],
            facecolor="None",
            edgecolor='r',
            s=12,
            linewidth=1,
            label='satellite galaxies'
            )
        ax.scatter(xg[0], yg[0],
            facecolor="None",
            edgecolor='lime',
            #  edgecolor='yellow',
            s=12,
            linewidth=1,
            label='central galaxy'
            )


        ax.legend(loc='lower right', framealpha=0.5, fontsize=15)







    #----------------------
    # Ax specific tweaks
    #----------------------


    title = "Projection of DM particles for halo "+str(halonr)+', z={0:5.3f}'.format(abs(1.0/aexp-1))
    ax1.set_title(title, fontsize=16)

    title = "Galaxies with host DM clumps"
    ax2.set_title(title, fontsize=16)

    title = "Galaxies with host DM clumps and orphan galaxies"
    ax3.set_title(title, fontsize=16)


    if horizontal:
        plt.subplots_adjust(left=0.01, right=0.99, top=0.97, bottom=0.01,wspace=0.01)
        #  plt.figtext(.02, .03, 'z={0:5.3f}'.format(z))
    else:
        plt.subplots_adjust(left=0.00, right=1,top=0.98,bottom=0.01,hspace=0.07)

    outfname = 'halostats_image_output_'+nr+'-halo-'+str(halonr)+'.pdf'
    plt.savefig(outfname, format='pdf')
    #  print("Saved fig", outfname)
    #  outfname = 'halostats_image_output_'+nr+'-halo-'+str(halonr)+'.png'
    #  plt.savefig(outfname, format='png')
    plt.close()




#==========================================
def plot_halo(dirnr, halonr):
#==========================================
    """
    Plots pdf image of halo and galaxies.
    Needs output_XXXXX number and halo number.
    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    nr = str(dirnr).zfill(5)

    #  print("Plotting ", halonr, "from output_"+nr)

    #---------------------------
    # get halo image data
    #---------------------------

    fname = "output_"+nr+'/density_image_new-halo-'+str(halonr)+'.dat'
    f = ff.FortranFile(fname)
    nc = np.asscalar(f.readInts())
    aexp   = np.asscalar(f.readReals('d'))
    unit_l = np.asscalar(f.readReals('d'))
    xmax   = np.asscalar(f.readReals('d'))
    xmin   = np.asscalar(f.readReals('d'))
    ymax   = np.asscalar(f.readReals('d'))
    ymin   = np.asscalar(f.readReals('d'))
    zmax   = np.asscalar(f.readReals('d'))
    zmin   = np.asscalar(f.readReals('d'))
    imageXY  = f.readReals('d')
    imageYZ  = f.readReals('d')
    imageXZ  = f.readReals('d')


    #---------------------------
    # get halo metadata
    #---------------------------

    def getfloat(f):
        line = f.readline()
        line = line.strip()
        return float(line)

    if use_partcenter:
        fname = "output_"+nr+'/radial-profile-partcenter-metadata-'+str(halonr)+'.dat'
    else:
        fname = "output_"+nr+'/radial-profile-metadata-'+str(halonr)+'.dat'
    f = open(fname)
    f.readline() # skip comment
    xc = getfloat(f)
    yc = getfloat(f)
    zc = getfloat(f)
    f.readline() # skip comment
    r200 = getfloat(f)
    f.readline() # skip comment
    rmax = getfloat(f)
    f.readline() # skip comment
    m200 = getfloat(f)
    f.readline() # skip comment
    unit_l = getfloat(f)
    f.readline() # skip comment
    unit_l_Mpc = getfloat(f)
    f.close()

    if not use_internal_units:
        xmin *= unit_l_Mpc
        xmax *= unit_l_Mpc
        ymin *= unit_l_Mpc
        ymax *= unit_l_Mpc
        zmin *= unit_l_Mpc
        zmax *= unit_l_Mpc


    #  print "r200: {0:12.3e}/{1:12.3e} Mpc".format(r200, r200*unit_l_Mpc)
    #  print "rmax: {0:12.3e}/{1:12.3e} Mpc".format(rmax, rmax*unit_l_Mpc)





    for image in [imageXY, imageYZ, imageXZ]:
        minval=np.min(image[image>0])
        maxval = np.max(image)
        minval/=10
        # cheat to not have white spaces with lognorm on image
        #  image = np.swapaxes(image, 0, 1)
        image[image==0] = minval
    imageXY = imageXY.reshape((nc,nc))
    imageYZ = imageYZ.reshape((nc,nc))
    imageXZ = imageXZ.reshape((nc,nc))


    plt.close("all") # safety measure

    # Plot image
    if horizontal:
        fig = plt.figure(figsize=(18,6.2))
        ax1 = fig.add_subplot(1,3,1, aspect="equal")
        ax2 = fig.add_subplot(1,3,2, aspect="equal")
        ax3 = fig.add_subplot(1,3,3, aspect="equal")
    else:
        fig = plt.figure(figsize=(6.2,18))
        ax1 = fig.add_subplot(3,1,1, aspect="equal")
        ax2 = fig.add_subplot(3,1,2, aspect="equal")
        ax3 = fig.add_subplot(3,1,3, aspect="equal")



    #------------------------------
    # Plot Halo Backgrounds
    #------------------------------

    def plot_data(ax, data, xmin, xmax, ymin, ymax):

        ax.imshow(
            data, 
            origin='lower',
            extent=(xmin, xmax, ymin, ymax),
            norm=colors.LogNorm(),
            interpolation='Kaiser',
            cmap='Blues',
            vmin=minval/100, vmax=maxval
        )

    def plot_circles(ax, xc, yc, r200, rmax):

        ax.scatter([xc], [yc], c='r', s=10)

        rrmax_pix = ax.transData.transform(np.vstack([rmax, rmax]).T) - ax.transData.transform(np.vstack([np.zeros(1), np.zeros(1)]).T)
        rmax_pix, _ = rrmax_pix.T

        size_pt_rmax = (2 * rmax_pix / fig.dpi * 72) ** 2
        s1 = ax.scatter([xc], [yc], facecolor='none', edgecolor='k', s=0)
        s1.set_sizes(size_pt_rmax)

        rr200_pix = ax.transData.transform(np.vstack([r200, r200]).T) - ax.transData.transform(np.vstack([np.zeros(1), np.zeros(1)]).T)
        r200_pix, _ = rr200_pix.T

        size_pt_r200 = (2 * r200_pix / fig.dpi * 72) ** 2
        s2 = ax.scatter([xc], [yc], facecolor='none', edgecolor='r', s=0)
        s2.set_sizes(size_pt_r200)


    plot_data(ax1, imageXY, xmin, xmax, ymin, ymax)
    plot_data(ax2, imageYZ, ymin, ymax, zmin, zmax)
    plot_data(ax3, imageXZ, xmin, xmax, zmin, zmax)

    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax2.set_xlim(ymin, ymax)
    ax2.set_ylim(zmin, zmax)
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(zmin, zmax)

    fig.canvas.draw()

    plot_circles(ax1, xc, yc, r200, rmax)
    plot_circles(ax2, yc, zc, r200, rmax)
    plot_circles(ax3, xc, zc, r200, rmax)


    #  for ax in [ax1, ax2, ax3]:
    for ax in [ax1]:

        #  ax.set_xticklabels([])
        #  ax.set_yticklabels([])
        #  ax.set_xticks([])
        #  ax.set_yticks([])

        dx = xmax - xmin
        xs = xmin + 0.05*dx
        xe = xmin + 0.25*dx
        ys = ymin + 0.05*dx
        ye = ymin + 0.07*dx
        xm = 0.5*(xs+xe)

        plt.sca(ax)
        plt.annotate("",
            (xs, ys), (xe, ys),
            arrowprops={
                'arrowstyle':'<->',
                'color':'m',
                'linewidth':'1'
                }
            )

        dx = xe - xs
        arrowlength = dx
        if use_internal_units:
            arrowlength *= unit_l_Mpc
            #  print "dx {0:12.3e} {1:12.3e}".format(dx, dx*unit_l_Mpc)
        plt.annotate("{0:3.1f} Mpc".format(arrowlength),
            xy=(xm, ye),
            xycoords='data',
            color='m',
            horizontalalignment='center',
            fontsize=15
            )







    #----------------------
    # Ax specific tweaks
    #----------------------


    #  title = "Projection of DM particles for halo "+str(halonr)+', z={0:5.3f}'.format(abs(1.0/aexp-1))
    #  ax1.set_title(title, fontsize=16)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")

    #  title = "Galaxies with host DM clumps"
    #  ax2.set_title(title, fontsize=16)
    ax2.set_xlabel("y")
    ax2.set_ylabel("z")

    #  title = "Galaxies with host DM clumps and orphan galaxies"
    #  ax3.set_title(title, fontsize=16)
    ax3.set_xlabel("x")
    ax3.set_ylabel("z")

    fig.suptitle("Halo "+str(halonr))


    plt.tight_layout()

    outfname = 'halo-density-image_output_'+nr+'-halo-'+str(halonr)+'.png'
    plt.savefig(outfname, format='png')
    #  print("Saved fig", outfname)
    plt.close()








#===================
def main():
#===================

    from sys import argv

    outputnr = argv[1]
    halonr = argv[2]


    if plot_galaxies:
        plot_halo_and_galaxies(outputnr, halonr)
    else:
        plot_halo(outputnr, halonr)





#==============================
if __name__ == "__main__":
#==============================
    main()

