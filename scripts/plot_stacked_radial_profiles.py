#!/usr/bin/env python3

#===================================================================
# stack radial profiles generated by plot_radial_profiles.py
# in combination with get_halostats.sh, plot them, and
# fit NFW profiles on top of it
#
# the script is looking for all available files named
# radial-profiles-data-XXXXX-halo-<int>.pkl
#
# This script needs one command line argument: the output number
# to work with
#===================================================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from scipy.optimize import curve_fit
import pickle
import os

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif',
    'serif':['Computer Modern Roman'],
    'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)
rc("figure", **{"dpi":200})



#==========================================
def sigma_nfw(r, c, rho0):
#==========================================
    """
    get NFW surface density profile
    assumes r is a numpy array
    following Bartelmann 1996, https://arxiv.org/pdf/astro-ph/9602053.pdf
    """
    Rs = 1./ c
    x = r / Rs

    f = np.zeros(x.shape)
    f[x > 1] = 1 - 2/np.sqrt(x[x>1]**2 - 1) * np.arctan(np.sqrt((x[x>1] - 1)/(x[x>1] + 1)))
    f[x < 1] = 1 - 2/np.sqrt(1 - x[x<1]**2) * np.arctanh(np.sqrt((1 - x[x<1])/(1 + x[x<1])))
    f[x == 1] = 0.

    return  2 * rho0 * Rs / (x**2 - 1) * f



#================================
def fit_nfw(xdata, ydata):
#================================
    """
    Get a fit to NFW surface density profile

    Returns: 
        (c, rho): Parameters to be fitted
    """
    
    cmax = 20.
    rhoguess = 0.5 * (ydata[0] + ydata[-1])
    opt, cov = curve_fit(sigma_nfw, xdata, ydata, p0 = [0.5, rhoguess], bounds=(0, [cmax, 1e30]), maxfev=10000000)
    c, rho = opt
    if abs(c/cmax - 1) < 0.01:
        # try again with reduced dataset
        if xdata.shape[0] > 6:
            print("restarting fit")
            return fit_nfw(xdata[1:-1], ydata[1:-1])
    return c, rho



#==========================================
def plot_stacked_radial_profile(outputnr):
#==========================================
    """
    Plots radial profile.
    """

    import matplotlib.pyplot as plt
    #  import matplotlib.colors as colors


    # find all ellegible pickle files
    dirlist = os.listdir(os.getcwd())
    filelist = []
    #  namestart = 'radial-profiles-data-'+outputnr.zfill(5)
    namestart = 'radial-profiles-data-for-stacked-profile-'+outputnr.zfill(5)
    for f in dirlist:
        if f.startswith(namestart) and f.endswith(".pkl"):
            filelist.append(f)

    # read thresholds
    f = open(filelist[0], 'rb')
    thresholds = pickle.load(f)
    threshold_labels = pickle.load(f)
    f.close()

    # open all files and skip to data
    filepointers = [open(fname, "rb") for fname in filelist]
    surface_densities_g = [] # list of lists of corresponding data. Each loaded file has number of thresholds elements in its list.
    surface_densities_o = []
    number_densities_g = []
    number_densities_o = []
    lin_surface_densities_g = []
    lin_surface_densities_o = []
    lin_number_densities_g = []
    lin_number_densities_o = []

    for p in filepointers:
        _ = pickle.load(p) # skip thresholds
        _ = pickle.load(p) # skip threshold labels
        surface_densities_g.append(pickle.load(p))
        surface_densities_o.append(pickle.load(p))
        number_densities_g.append(pickle.load(p))
        number_densities_o.append(pickle.load(p))
        bin_centres = pickle.load(p) 
        lin_surface_densities_g.append(pickle.load(p))
        lin_surface_densities_o.append(pickle.load(p))
        lin_number_densities_g.append(pickle.load(p))
        lin_number_densities_o.append(pickle.load(p))
        lin_bin_number_centres = pickle.load(p)
        lin_bin_mass_centres = pickle.load(p)


    # stack data together
    stacked_surface_densities_g = []
    stacked_surface_densities_o = []
    stacked_number_densities_g = []
    stacked_number_densities_o = []

    stacked_lin_surface_densities_g = []
    stacked_lin_surface_densities_o = []
    stacked_lin_number_densities_g = []
    stacked_lin_number_densities_o = []

    stacked_surface_density_fit_g = []
    stacked_surface_density_fit_o = []
    stacked_number_density_fit_g = []
    stacked_number_density_fit_o = []

    for i, t in enumerate(thresholds):
        sg = None
        for halo in surface_densities_g:
            data = halo[i]
            if sg is None:
                sg = np.copy(data)
            else:
                sg += data
        stacked_surface_densities_g.append(sg)

    for i, t in enumerate(thresholds):
        sg = None
        for halo in lin_surface_densities_g:
            data = halo[i]
            if sg is None:
                sg = np.copy(data)
            else:
                sg += data
        stacked_lin_surface_densities_g.append(sg)

        mask = sg > 0
        ydata = sg[mask]
        xdata = lin_bin_mass_centres[mask]
        stacked_surface_density_fit_g.append(fit_nfw(xdata, ydata))

    for i, t in enumerate(thresholds):
        so = None
        for halo in surface_densities_o:
            data = halo[i]
            if so is None:
                so = np.copy(data)
            else:
                so += data
        stacked_surface_densities_o.append(so)

    for i, t in enumerate(thresholds):
        so = None
        for halo in lin_surface_densities_o:
            data = halo[i]
            if so is None:
                so = np.copy(data)
            else:
                so += data
        stacked_lin_surface_densities_o.append(so)

        mask = so > 0
        ydata = so[mask]
        xdata = lin_bin_mass_centres[mask]
        stacked_surface_density_fit_o.append(fit_nfw(xdata, ydata))


    for i, t in enumerate(thresholds):
        ng = None
        for halo in number_densities_g:
            data = halo[i]
            if ng is None:
                ng = np.copy(data)
            else:
                ng += data
        stacked_number_densities_g.append(ng)

    for i, t in enumerate(thresholds):
        ng = None
        for halo in lin_number_densities_g:
            data = halo[i]
            if ng is None:
                ng = np.copy(data)
            else:
                ng += data
        stacked_lin_number_densities_g.append(ng)

        mask = ng > 0
        ydata = ng[mask]
        xdata = lin_bin_number_centres[mask]
        stacked_number_density_fit_g.append(fit_nfw(xdata, ydata))


    for i, t in enumerate(thresholds):
        no = None
        for halo in number_densities_o:
            data = halo[i]
            if no is None:
                no = np.copy(data)
            else:
                no += data
        stacked_number_densities_o.append(no)

    for i, t in enumerate(thresholds):
        no = None
        for halo in lin_number_densities_o:
            data = halo[i]
            if no is None:
                no = np.copy(data)
            else:
                no += data
        stacked_lin_number_densities_o.append(no)

        mask = no > 0
        ydata = no[mask]
        xdata = lin_bin_number_centres[mask]
        stacked_number_density_fit_o.append(fit_nfw(xdata, ydata))



    fig = plt.figure(1, figsize=(10, 10))
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']
    alpha = 0.8
    fitalpha = 0.8
    nonstackedalpha = 0.3
    nonstackedcolor = 'grey'
    nonstackedfactor = 1e-3


    # plot number densities
    ax1 = fig.add_subplot(221)
    ax1.set_title("Surface Number Densities of Galaxies with Hosts")
    ax1.set_ylabel(r"$\Sigma_n$ [$r_{200}\ ^{-2}$]")

    for i, t in enumerate(thresholds):
        number_density_g = stacked_number_densities_g[i]
        mask = number_density_g >= 0
        ax1.loglog(
            bin_centres[mask], 
            number_density_g[mask], 
            label="galaxies with hosts "+threshold_labels[i], 
            c = colors[i], 
            alpha = alpha, 
            )

        c, rho = stacked_number_density_fit_g[i]
        ax1.loglog(
            bin_centres,
            sigma_nfw(bin_centres, c, rho),
            "-.",
            c=colors[i],
            label="fit with c={0:.2f}".format(c),
            alpha = fitalpha,
            )

        #  for nd in number_densities_g:
        #      ax1.loglog(
        #                  bin_centres,
        #                  nd[0] * nonstackedfactor,
        #                  c = nonstackedcolor,
        #                  alpha = nonstackedalpha,
        #              )
        #

    # make an emptly plot with the legend

    # plot number densities
    ax2 = fig.add_subplot(222)
    ax2.set_title("Surface Number Densities of Galaxies with Hosts and Orphans")
    ax2.set_ylabel(r"$\Sigma_n$ [$r_{200}\ ^{-2}$]")

    for i, t in enumerate(thresholds):
        number_density_o = stacked_number_densities_o[i]
        mask = number_density_o >= 0
        ax2.loglog(
                    bin_centres[mask], 
                    number_density_o[mask], 
                    label="all galaxies, including orphans "+threshold_labels[i],
                    c = colors[i],
                    alpha = alpha, 
                    )

        c, rho = stacked_number_density_fit_o[i]
        ax2.loglog(
            bin_centres,
            sigma_nfw(bin_centres, c, rho),
            "-.",
            c=colors[i],
            label="fit with c={0:.2f}".format(c),
            alpha = fitalpha,
            )

        #  for nd in number_densities_o:
        #      ax2.loglog(
        #                  bin_centres,
        #                  nd[0] * nonstackedfactor,
        #                  c = nonstackedcolor,
        #                  alpha = nonstackedalpha,
        #              )
        #

    # make an emptly plot with the legend
    ax3 = fig.add_subplot(223)
    ax3.set_title("Surface Mass Densities of Galaxies with Hosts")
    ax3.set_ylabel(r"$\Sigma_m$ [$M_\odot$ $r_{200}\ ^{-2}$]")
    for i, t in enumerate(thresholds):
        surface_density_g = stacked_surface_densities_g[i]
        mask = surface_density_g > 0
        ax3.loglog(
                    bin_centres[mask], 
                    surface_density_g[mask], 
                    label="galaxies with hosts "+threshold_labels[i], 
                    c = colors[i], 
                    alpha = alpha, 
                    )

        c, rho = stacked_surface_density_fit_g[i]
        ax3.loglog(
            bin_centres,
            sigma_nfw(bin_centres, c, rho),
            "-.",
            c=colors[i],
            label="fit with c={0:.2f}".format(c),
            alpha = fitalpha,
            )

        #  for sd in surface_densities_g:
        #      ax3.loglog(
        #                  bin_centres,
        #                  sd[0] * nonstackedfactor,
        #                  c = nonstackedcolor,
        #                  alpha = nonstackedalpha,
        #              )
        #

    # make an emptly plot with the legend


    ax4 = fig.add_subplot(224)
    ax4.set_title("Surface Mass Densities of Galaxies with Hosts and orphans")
    ax4.set_ylabel(r"$\Sigma_m$ [$M_\odot$ $r_{200}\ ^{-2}$]")
    for i, t in enumerate(thresholds):
        surface_density_o = stacked_surface_densities_o[i]
        mask = surface_density_o > 0

        ax4.loglog(
            bin_centres[mask], 
            surface_density_o[mask], 
            label="including orphans "+threshold_labels[i], 
            c = colors[i], 
            alpha = alpha, 
            )

        c, rho = stacked_surface_density_fit_o[i]
        ax4.loglog(
            bin_centres,
            sigma_nfw(bin_centres, c, rho),
            "-.",
            c=colors[i],
            label="fit with c={0:.2f}".format(c),
            alpha = fitalpha,
            )

        #  for sd in surface_densities_o:
        #      ax4.loglog(
        #                  bin_centres,
        #                  sd[0] * nonstackedfactor,
        #                  c = nonstackedcolor,
        #                  alpha = nonstackedalpha,
        #              )
        #


    for ax in fig.axes:
        ax.legend(loc="lower left")
        ax.set_xlabel(r"$r/r_{200}$")
        ax.grid()



    # save image
    fname = "radial-profiles-stacked-"+outputnr+".png"
    plt.tight_layout()
    #  plt.show()
    plt.savefig(fname)
    print("saved figure", fname)
    plt.close()



    return






#===================
def main():
#===================

    from sys import argv

    outputnr = argv[1]

    plot_stacked_radial_profile(outputnr)





#==============================
if __name__ == "__main__":
#==============================
    main()

