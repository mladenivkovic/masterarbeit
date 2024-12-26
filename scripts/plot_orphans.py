#!/usr/bin/python2
# need python2 here for fortranfile module

#===========================================
# Plots halos from output made by
# get_orphan_image_data.f03, which is compiled and
# called by get_orphan_plots.sh
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
mthresh = 10.



#==========================================
def create_plots(dirnr, halonr):
#==========================================
    """
    Plots pdf image of halo and galaxies.
    Needs output_XXXXX number and halo number.
    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    nr = str(dirnr).zfill(5)

    #---------------------------
    # get halo image data
    #---------------------------

    fname_in = "orphan_plots/orphan_plots-"+nr+"/orphan_plot-halo-"+str(halonr)+".dat"
    fname_out = fname_in[:-3] + "png"

    f = ff.FortranFile(fname_in)
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
    
    ng = np.asscalar(f.readInts())
    xg = f.readReals('d')
    yg = f.readReals('d')
    zg = f.readReals('d')
    mg = f.readReals('d')
    logmg = np.log10(mg)

    no = np.asscalar(f.readInts())
    xo = f.readReals('d')
    yo = f.readReals('d')
    zo = f.readReals('d')
    mo = f.readReals('d')
    logmo = np.log10(mo)

    #  massmin = min(logmg.min(), logmo.min())
    #  massmax = min(logmg.max(), logmo.max())
    massmin = 9.8
    massmax = 13.0

    # select galaxies above threshold
    mask = logmg > mthresh
    logmg = logmg[mask]
    xg = xg[mask]
    yg = yg[mask]
    zg = zg[mask]
    ng = xg.shape[0]

    mask = logmo > mthresh
    logmo = logmo[mask]
    xo = xo[mask]
    yo = yo[mask]
    zo = zo[mask]
    no = xo.shape[0]


    for image in [imageXY, imageYZ, imageXZ]:
        minval=np.min(image[image>0])
        maxval = np.max(image)
    minval/=10

    for image in [imageXY, imageYZ, imageXZ]:
        # cheat to not have white spaces with lognorm on image
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
        fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)
    else:
        fig = plt.figure(figsize=(6.2,18))
        ax1 = fig.add_subplot(3,1,1, aspect="equal")
        ax2 = fig.add_subplot(3,1,2, aspect="equal")
        ax3 = fig.add_subplot(3,1,3, aspect="equal")
        fig.subplots_adjust(left=0.02, bottom=0.06, right=0.95, top=0.94, wspace=0.05)



    #------------------------------
    # Plot Halo Backgrounds
    #------------------------------

    def plot_background(ax, data, xmin, xmax, ymin, ymax):

        ax.imshow(
            data, 
            origin='lower',
            extent=(xmin, xmax, ymin, ymax),
            norm=colors.LogNorm(),
            interpolation='Kaiser',
            cmap='Blues',
            vmin=minval/10, vmax=maxval
        )

    def plot_galaxies(fig, ax, x, y, m):

        sc = ax.scatter(
            x, y, c=m, 
            cmap='afmhot',
            vmin=massmin, 
            vmax=massmax, 
            s = 120,
            marker = '*', 
            lw = 0.2, 
            alpha = 1., 
            )

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(sc, cax=cax)

        return

    def plot_orphans(fig, ax, x, y, m):

        sc = ax.scatter(
            x, y, c=m, 
            cmap='afmhot',
            vmin=massmin, 
            vmax=massmax, 
            s = 20,
            marker = 'o', 
            lw = 0.2, 
            alpha = 1., 
            )

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(sc, cax=cax)

        return

    plot_background(ax1, imageXY, xmin, xmax, ymin, ymax)
    plot_galaxies(fig, ax1, xg, yg, logmg)
    plot_orphans(fig, ax1, xo, yo, logmo)
    plot_background(ax2, imageYZ, ymin, ymax, zmin, zmax)
    plot_galaxies(fig, ax2, yg, zg, logmg)
    plot_orphans(fig, ax2, yo, zo, logmo)
    plot_background(ax3, imageXZ, xmin, xmax, zmin, zmax)
    plot_galaxies(fig, ax3, xg, zg, logmg)
    plot_orphans(fig, ax3, xo, zo, logmo)



    ax1.set_xlim(xmin, xmax)
    ax1.set_ylim(ymin, ymax)
    ax2.set_xlim(ymin, ymax)
    ax2.set_ylim(zmin, zmax)
    ax3.set_xlim(xmin, xmax)
    ax3.set_ylim(zmin, zmax)





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

    fig.suptitle(
        r"log($M_*$) for galaxies (stars; " + \
        str(ng) + ") and orphans (dots; " + str(no) + \
        ") with log($M_*$) $>$ {0:.1f}".format(mthresh)+" for halo " + str(halonr))


    plt.tight_layout()

    plt.savefig(fname_out, format='png')
    #  print("Saved fig", outfname)
    plt.close()








#===================
def main():
#===================

    from sys import argv

    outputnr = argv[1]
    halonr = argv[2]

    create_plots(outputnr, halonr)





#==============================
if __name__ == "__main__":
#==============================
    main()

