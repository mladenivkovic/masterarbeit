#!/usr/bin/env python3

#===========================================
# Plot radial profiles and fit concentration
# c to NFW surface density profiles for 
# galaxies including or excluding orphans
#
# Made to work with get_halostats.sh output
# which in return calls get_halostats.f03
#
# needs halo number and output number as 
# cmdline arguments. output number needs to
# be zero paddded:
# plot_radial_profiles.py 0067 123456
#===========================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from scipy.optimize import curve_fit
import pickle

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif',
    'serif':['Computer Modern Roman'],
    'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)
rc("figure", **{"dpi":200})

print_table_data = False
table_data_r200_cutoff = True
make_plots = True # skip actual plotting if you only want output for stacked profiles


nlogbins = 20
nlinbins = 200
rgal = 10*1e-3 # 10 kpc

alpha = 0.8
fitalpha = 0.8

Mpc = 3.086e24


#==========================================
def sigma_nfw(r, c, rho0):
#==========================================
    """
    get NFW surface density profile
    assumes r is a numpy array
    following Bartelmann 1996, https://arxiv.org/pdf/astro-ph/9602053.pdf
    """
    Rs = 1. / c  # divide r200 by itself because you want everything in units of r200
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
            #  print("restarting fit")
            return fit_nfw(xdata[1:-1], ydata[1:-1])
    return c, rho



#======================================================================
def plot_radial_profile(outputnr, halonr):
#======================================================================
    """
    Plots radial profile.
    Create a simple radial histogram.
    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    nr = str(outputnr).zfill(5)

    def getfloat(f):
        line = f.readline()
        line = line.strip()
        return float(line)



    #  print("Computing radial profile")

    #---------------------------
    # get halo metadata
    #---------------------------

    fname = "output_"+nr+'/radial-profile-partcenter-metadata-'+str(halonr)+'.dat'
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

    r200 *= unit_l_Mpc
    rmax *= unit_l_Mpc


    #-------------------------------
    # Read in data
    #-------------------------------

    # get galaxy data
    fname = "output_"+nr+'/profile-galaxies-'+str(halonr)+'.dat'
    mg, xg, yg, zg = np.loadtxt(fname, dtype='float', skiprows=1, unpack=True)

    # get orphans data
    fname = "output_"+nr+'/profile-orphans-'+str(halonr)+'.dat'
    mo, xo, yo, zo = np.loadtxt(fname, dtype='float', skiprows=1, unpack=True)

    #  print("m200: {0:12.3e}, central mass: {1:12.3e}, total stellar mass: {2:12.3e}".format(m200, mg[0], mg[1:].sum()))

    xc *= unit_l_Mpc
    yc *= unit_l_Mpc
    zc *= unit_l_Mpc


    #-------------------------------------
    # Prepare data for furher treatment
    #-------------------------------------

    # remove central galaxy from list
    mcentral = mg[0]
    mg = mg[1:]
    xg = xg[1:]
    yg = yg[1:]
    zg = zg[1:]

    # get radii from central galaxy
    # ignore z coordinate: assume we're on the x-y plane
    rg = (xg - xc)**2 + (yg - yc)**2 #+ (zg - zc)**2
    ro = (xo - xc)**2 + (yo - yc)**2 #+ (zo - zc)**2
    # TODO: periodicity?
    if ((xg - xc) > 50.).any():
        print("periodicity error?")
        print(xg)
        quit()
    if ((xg - xc) < -50.).any():
        print("periodicity error?")
        quit()
    if ((yg - yc) > 50.).any():
        print("periodicity error?")
        quit()
    if ((yg - yc) < -50.).any():
        print("periodicity error?")
        quit()
    if ((zg - zc) > 50.).any():
        print("periodicity error?")
        quit()
    if ((zg - zc) < -50.).any():
        print("periodicity error?")
        quit()
    #  print("Max diff center:", (xg - xc).max(), (yg - yc).max(), (zg - zc).max())
    #  print("Max diff center orphans:", (xo - xc).max(), (yo - yc).max(), (zo - zc).max())

    # do everything in units of r200
    rg = np.sqrt(rg) / r200
    ro = np.sqrt(ro) / r200


    #  rmin = min(ro.min(), rg.min())
    #  rmin = max(rmin, 0.01)# don't go below 0.01 r_200
    #  rmax = max(ro.max(), rg.max())

    # set limits by hand
    rmin = 0.05
    rmax = 1.

    #  print("rmin in kpc:", rmin * r200 * 1000)

    #  print("r_min:", rmin, "r200")
    #  print("r_min:", rmin*r200*1000, "kpc") # default unit is Mpc

    logbins = np.logspace(np.log10(rmin), np.log10(rmax*2), num=nlogbins)

    #  linbins = np.linspace(rmin*0.9, rmax*1.1, num=nlinbins)
    # use varying bin widths to give more emphasis on outer regions of a halo
    # number density seems fine. Look to possibly improve mass density!
    linbins_number_density = np.hstack((np.linspace(rmin*0.99, 0.1, 10), np.linspace(0.11, rmax*1.01, 100)))
    linbins_mass_density = np.hstack((np.linspace(rmin*0.99, 0.1, 10), np.linspace(0.11, rmax*1.01, 100)))

    surface_densities_g = []
    surface_densities_o = []
    number_densities_g = []
    number_densities_o = []
    lin_surface_densities_g = []
    lin_surface_densities_o = []
    lin_number_densities_g = []
    lin_number_densities_o = []
    surface_density_fit_g = []
    surface_density_fit_o = []
    number_density_fit_g = []
    number_density_fit_o = []

    thresholds = [0, 1e9, 1e10]
    threshold_labels = ["no threshold", r"$M > 10^9 M_\odot$", r"$M > 10^{10} M_\odot$"]


    # open dump file
    fname = "radial-profiles-data-for-stacked-profile-"+outputnr.zfill(5)+"-halo-"+halonr+".pkl"
    dump = open(fname, 'wb')

    # dump thresholds and labels
    pickle.dump(thresholds, dump)
    pickle.dump(threshold_labels, dump)

    if print_table_data:
        
        if table_data_r200_cutoff:
            mask_gals = rg <= 1.
            mask_orph = ro <= 1.
            ngals = mg[mask_gals].shape[0]
            norph = mo[mask_orph].shape[0]
            ngals1 = np.count_nonzero(mg[mask_gals] > thresholds[1])
            ngals2 = np.count_nonzero(mg[mask_gals] > thresholds[2])
            norph1 = np.count_nonzero(mo[mask_orph] > thresholds[1])
            norph2 = np.count_nonzero(mo[mask_orph] > thresholds[2])
        else:
            ngals = mg.shape[0]
            norph = mo.shape[0]
            ngals1 = np.count_nonzero(mg > thresholds[1])
            ngals2 = np.count_nonzero(mg > thresholds[2])
            norph1 = np.count_nonzero(mo > thresholds[1])
            norph2 = np.count_nonzero(mo > thresholds[2])
        print("{0:8d} & {1:12.2e} & {2:12.2e} & {3:12.2e} & Count            & {4:12d} & {5:12d} & {6:12d} & {7:12d} & {8:12d} & {9:12d}\\\\".format(
                int(halonr), m200, r200, mcentral, ngals, norph, ngals1, norph1, ngals2, norph2))

        if table_data_r200_cutoff:
            mask_gals = rg <= 1.
            mask_orph = ro <= 1.
            mgals = mg[mask_gals].sum()
            morph = mo[mask_orph].sum()
            mgals1 = np.sum(mg[mask_gals][mg[mask_gals] > thresholds[1]])
            mgals2 = np.sum(mg[mask_gals][mg[mask_gals] > thresholds[2]])
            morph1 = np.sum(mo[mask_orph][mo[mask_orph] > thresholds[1]])
            morph2 = np.sum(mo[mask_orph][mo[mask_orph] > thresholds[2]])
        else:
            mgals = mg.sum()
            morph = mo.sum()
            mgals1 = np.sum(mg[mg > thresholds[1]])
            mgals2 = np.sum(mg[mg > thresholds[2]])
            morph1 = np.sum(mo[mo > thresholds[1]])
            morph2 = np.sum(mo[mo > thresholds[2]])
        print("{0:8} & {0:12} & {0:12} & {0:12} & Mass [$M_\\odot$] & {1:12.2e} & {2:12.2e} & {3:12.2e} & {4:12.2e} & {5:12.2e} & {6:12.2e}\\\\".format(
                " ", mgals, morph, mgals1, morph1, mgals2, morph2))
        print("\\hline")



    for t in thresholds:
        #------------------------
        # Histogram Data
        #------------------------

        choice_g = mg > t
        choice_o = mo > t
        #  print("Threshold", t, "keeping",
        #          np.count_nonzero(choice_g), "/", choice_g.shape[0], "galaxies;",
        #          np.count_nonzero(choice_o), "/", choice_o.shape[0], "orphans;")

        number_g, bin_edges = np.histogram(rg[choice_g], bins=logbins)
        mass_g, bin_edges = np.histogram(rg[choice_g], bins=logbins, weights=mg[choice_g])
        number_o, bin_edges = np.histogram(ro[choice_o], bins=logbins, )
        mass_o, bin_edges = np.histogram(ro[choice_o], bins=logbins, weights=mo[choice_o])

        # add galaxies to orphans
        mass_o += mass_g
        number_o += number_g
        bin_centres = (bin_edges[:-1] + bin_edges[1:]) * 0.5

        surfaces = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)
        surface_density_g = mass_g / surfaces
        number_density_g = number_g / surfaces
        surface_density_o = mass_o / surfaces
        number_density_o = number_o / surfaces

        surface_densities_g.append(surface_density_g)
        surface_densities_o.append(surface_density_o)
        number_densities_g.append(number_density_g)
        number_densities_o.append(number_density_o)


        #------------------------------------
        # Histogram data lineraly as well
        #------------------------------------

        lin_number_g, lin_bin_number_edges = np.histogram(rg[choice_g], bins=linbins_number_density)
        lin_mass_g, lin_bin_mass_edges = np.histogram(rg[choice_g], bins=linbins_mass_density, weights=mg[choice_g])
        lin_number_o, lin_bin_number_edges = np.histogram(ro[choice_o], bins=linbins_number_density, )
        lin_mass_o, lin_bin_mass_edges = np.histogram(ro[choice_o], bins=linbins_mass_density, weights=mo[choice_o])

        # add galaxies to orphans
        lin_mass_o += lin_mass_g
        lin_number_o += lin_number_g
        lin_bin_number_centres = (lin_bin_number_edges[:-1] + lin_bin_number_edges[1:]) * 0.5
        lin_bin_mass_centres = (lin_bin_mass_edges[:-1] + lin_bin_mass_edges[1:]) * 0.5

        lin_number_surfaces = np.pi * (lin_bin_number_edges[1:]**2 - lin_bin_number_edges[:-1]**2)
        lin_mass_surfaces = np.pi * (lin_bin_mass_edges[1:]**2 - lin_bin_mass_edges[:-1]**2)
        lin_surface_density_g = lin_mass_g / lin_mass_surfaces
        lin_number_density_g = lin_number_g / lin_number_surfaces
        lin_surface_density_o = lin_mass_o / lin_mass_surfaces
        lin_number_density_o = lin_number_o / lin_number_surfaces

        lin_surface_densities_g.append(lin_surface_density_g)
        lin_surface_densities_o.append(lin_surface_density_o)
        lin_number_densities_g.append(lin_number_density_g)
        lin_number_densities_o.append(lin_number_density_o)



        #------------------------
        # Get Fits
        #------------------------
        
        if make_plots:
            mask = lin_number_density_o > 0
            ydata = lin_number_density_o[mask]
            xdata = lin_bin_number_centres[mask]
            number_density_fit_o.append(fit_nfw(xdata, ydata))

            mask = lin_surface_density_o > 0
            ydata = lin_surface_density_o[mask]
            xdata = lin_bin_mass_centres[mask]
            surface_density_fit_o.append(fit_nfw(xdata, ydata))

            mask = lin_number_density_g > 0
            ydata = lin_number_density_g[mask]
            xdata = lin_bin_number_centres[mask]
            number_density_fit_g.append(fit_nfw(xdata, ydata))

            mask = lin_surface_density_g > 0
            ydata = lin_surface_density_g[mask]
            xdata = lin_bin_mass_centres[mask]
            surface_density_fit_g.append(fit_nfw(xdata, ydata))


    # dump actual data
    pickle.dump(surface_densities_g, dump)
    pickle.dump(surface_densities_o, dump)
    pickle.dump(number_densities_g, dump)
    pickle.dump(number_densities_o, dump)
    pickle.dump(bin_centres, dump)
    pickle.dump(lin_surface_densities_g, dump)
    pickle.dump(lin_surface_densities_o, dump)
    pickle.dump(lin_number_densities_g, dump)
    pickle.dump(lin_number_densities_o, dump)
    pickle.dump(lin_bin_number_centres, dump)
    pickle.dump(lin_bin_mass_centres, dump)
    dump.close()








    #----------------------
    # create plots
    #----------------------

    if make_plots:

        fig = plt.figure(figsize=(8,8))

        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']


        # plot number densities
        ax1 = fig.add_subplot(221)
        ax1.set_title("Surface Number Densities of Galaxies with Hosts")
        ax1.set_ylabel(r"$\Sigma_n$ [$r_{200}\ ^{-2}$]")

        for i, t in enumerate(thresholds):
            number_density_g = number_densities_g[i]
            mask = number_density_g >= 0
            ax1.loglog(
                bin_centres[mask], 
                number_density_g[mask], 
                label="galaxies with hosts "+threshold_labels[i], 
                c = colors[i], 
                alpha = alpha, 
                )

            c, rho = number_density_fit_g[i]
            ax1.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                )

        # plot number densities
        ax2 = fig.add_subplot(222)
        ax2.set_title("Surface Number Densities of Galaxies with Hosts and Orphans")
        ax2.set_ylabel(r"$\Sigma_n$ [$r_{200}\ ^{-2}$]")

        for i, t in enumerate(thresholds):
            number_density_o = number_densities_o[i]
            mask = number_density_o >= 0
            ax2.loglog(
                        bin_centres[mask], 
                        number_density_o[mask], 
                        label="all galaxies, including orphans "+threshold_labels[i],
                        c = colors[i],
                        alpha = alpha, 
                        )

            c, rho = number_density_fit_o[i]
            ax2.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                )


        # make an emptly plot with the legend
        ax3 = fig.add_subplot(223)
        ax3.set_title("Surface Mass Densities of Galaxies with Hosts")
        ax3.set_ylabel(r"$\Sigma_m$ [$M_\odot$ $r_{200}\ ^{-2}$]")
        for i, t in enumerate(thresholds):
            surface_density_g = surface_densities_g[i]
            mask = surface_density_g > 0
            ax3.loglog(
                        bin_centres[mask], 
                        surface_density_g[mask], 
                        label="galaxies with hosts "+threshold_labels[i], 
                        c = colors[i], 
                        alpha = alpha, 
                        )

            c, rho = surface_density_fit_g[i]
            ax3.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                )


        ax4 = fig.add_subplot(224)
        ax4.set_title("Surface Mass Densities of Galaxies with Hosts and orphans")
        ax4.set_ylabel(r"$\Sigma_m$ [$M_\odot$ $r_{200}\ ^{-2}$]")
        for i, t in enumerate(thresholds):
            surface_density_o = surface_densities_o[i]
            mask = surface_density_o > 0

            ax4.loglog(
                bin_centres[mask], 
                surface_density_o[mask], 
                label="including orphans "+threshold_labels[i], 
                c = colors[i], 
                alpha = alpha, 
                )

            c, rho = surface_density_fit_o[i]
            ax4.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                )


        for ax in fig.axes:
            ax.legend()
            ax.set_xlabel(r"$r/r_{200}$")
            ax.grid()




        # save image
        fname = "radial-profiles-"+outputnr.zfill(5)+"-halo-"+halonr+".png"
        plt.tight_layout()
        plt.savefig(fname)
        #  print("saved figure", fname)
        #  plt.show()
        plt.close()




    return





#===================
def main():
#===================

    from sys import argv

    outputnr = argv[1]
    halonr = argv[2]

    plot_radial_profile(outputnr, halonr)

    return





#==============================
if __name__ == "__main__":
#==============================
    main()


def plot_radial_profile_old_v2(outputnr, halonr):
    """
    Plots radial profile.
    """

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    nr = str(outputnr).zfill(5)

    def getfloat(f):
        line = f.readline()
        line = line.strip()
        return float(line)



    print("Computing radial profile")

    # get halo metadata

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
    f.close()



    # get galaxy data
    fname = "output_"+nr+'/profile-galaxies-'+str(halonr)+'.dat'
    mg, xg, yg, zg = np.loadtxt(fname, dtype='float', skiprows=1, unpack=True)

    # get orphans data
    fname = "output_"+nr+'/profile-orphans-'+str(halonr)+'.dat'
    mo, xo, yo, zo = np.loadtxt(fname, dtype='float', skiprows=1, unpack=True)


    print("Total stellar mass v1", (np.sum(mg) + np.sum(mo))*1e-12)
    print("Total stellar mass v2", (np.sum(mg[1:]) + np.sum(mo))*1e-12)

    #  print("setting central galaxy as center")
    #  xc = xg[0]
    #  yc = yg[0]
    #  zc = zg[0]

    # remove central galaxy from list
    mg = mg[1:]
    xg = xg[1:]
    yg = yg[1:]
    zg = zg[1:]


    # determine boundaries such that all projections have equal size
    xmin = min(xo.min(), xg.min())
    xmax = max(xo.max(), xg.max())
    ymin = min(yo.min(), yg.min())
    ymax = max(yo.max(), yg.max())
    zmin = min(zo.min(), zg.min())
    zmax = max(zo.max(), zg.max())

    #  print "x:", xmin, xmax
    #  print "y:", ymin, ymax
    #  print "z:", zmin, zmax

    dmax = max( abs(zmax - zc), 
                abs(zc - zmin),
                abs(ymax - yc), 
                abs(yc - ymin), 
                abs(xmax - xc),
                abs(xc - xmin))*1.01 # add a little extra for good taste

    xlo = xc - dmax
    ylo = yc - dmax
    zlo = zc - dmax
    xhi = xc + dmax
    yhi = yc + dmax
    zhi = zc + dmax

    dx = 2*dmax/nbins # dmax is only half box size
    dxhalf =  dx
    rmax_cells = (nbins + 0.5) * dx



    thresholds = [0, 1e9, 1e10]
    threshold_labels = ["no threshold", r"$M > 10^9 M_\odot$", r"$M > 10^{10} M_\odot$"]

    norphs_by_thresh = [None, None, None]
    ngals_by_thresh = [None, None, None]
    norph_fits = [None, None, None]
    ngals_fits = [None, None, None]

    morphs_by_thresh = [None, None, None]
    mgals_by_thresh = [None, None, None]
    morph_fits = [None, None, None]
    mgals_fits = [None, None, None]

    xg_all = xg
    yg_all = yg
    zg_all = zg
    mg_all = mg

    xo_all = xo
    yo_all = yo
    zo_all = zo
    mo_all = mo



    # open dump file
    fname = "radial-profiles-data-"+outputnr.zfill(5)+"-halo-"+halonr+".pkl"
    dump = open(fname, 'wb')

    # dump thresholds and labels
    pickle.dump(thresholds, dump)
    pickle.dump(threshold_labels, dump)



    def get_volume_fractions(ilo, ihi, x):
        """
        Get volume fractions inside cells with index ilo, ilo+1, ..., ihi
        ihi must be cell with the highest index where the galaxy is still overlapping
        x: galaxy position, where edge of box is origin, i.e. xg[i] - xlo

        returns: [(ilo, volume portion in this index), ..., (ihi, volume portion)]
        """
        dxi = []
        dxtot = 0
        if ilo == ihi:
            return [(ilo, 2*rgal)]

        for i in range(ilo, ihi+1):
            dx_lo = x - i*dx
            dx_hi =  x - (i+1)*dx
            if dx_lo > 0: # i*dx is below x
                if dx_lo > rgal:
                    dx_lo = rgal
            else:
                if abs(dx_lo) > rgal:
                    # stop if the lower edge of cell is already out of reach
                    break

            if dx_hi > 0: # (i+1)*dx is below x
                if dx_hi > rgal:
                    # upper boundary is above reach already, skip
                    continue
            else:
                if -dx_hi > rgal:
                    dx_hi = -rgal
                    
            # exception handling: ilo can be < 0, ihi may be > nbins
            current = i
            if current < 0:
                current = -i
            if current >= nbins:
                break

            dxi.append((current, abs(dx_hi - dx_lo)))
            dxtot += abs(dx_hi - dx_lo)

        # check that everything is in order
        if abs(dxtot / rgal) - 2  > 1e-5:
            print("dxtot / rgal =", dxtot/rgal, ", should be = 2. Something is wrong")
            quit()
        return dxi




    # Get number densities depending on Mass thresholds

    for I, mthresh in enumerate(thresholds):

        mask = mg_all > mthresh
        xg = xg_all[mask] - xlo
        yg = yg_all[mask] - ylo
        zg = zg_all[mask] - zlo
        mg = mg_all[mask]

        mask = mo_all > mthresh
        xo = xo_all[mask] - xlo
        yo = yo_all[mask] - ylo
        zo = zo_all[mask] - zlo
        mo = mo_all[mask]

        # create galaxy number projections

        ngals_xy,  xedges, yedges = np.histogram2d(xg, yg, bins=nbins, weights=None, range=((0, 2*dmax), (0, 2*dmax)))
        norphs_xy, xedges, yedges = np.histogram2d(xo, yo, bins=nbins, weights=None, range=((0, 2*dmax), (0, 2*dmax)))
        norphs_xy += ngals_xy

        ngals_yz,  xedges, yedges = np.histogram2d(yg, zg, bins=nbins, weights=None, range=((0, 2*dmax), (0, 2*dmax)))
        norphs_yz, xedges, yedges = np.histogram2d(yo, zo, bins=nbins, weights=None, range=((0, 2*dmax), (0, 2*dmax)))
        norphs_yz += ngals_yz

        ngals_zx,  xedges, yedges = np.histogram2d(zg, xg, bins=nbins, weights=None, range=((0, 2*dmax), (0, 2*dmax)))
        norphs_zx, xedges, yedges = np.histogram2d(zo, xo, bins=nbins, weights=None, range=((0, 2*dmax), (0, 2*dmax)))
        norphs_zx += ngals_zx



        # create galaxy mass projections
        mgals_xy = np.zeros((nbins, nbins), dtype='float')
        mgals_yz = np.zeros((nbins, nbins), dtype='float')
        mgals_zx = np.zeros((nbins, nbins), dtype='float')

        area = dx**2
        print("dx = ", 1000*dx, "kpc")
        for g in range(mg.shape[0]):

            # NGP interpolation
            #----------------------------------
            #  i = int(xg[g] / dx)
            #  j = int(yg[g] / dx)
            #  k = int(zg[g] / dx)
            #  mgals_xy[i, j] += mg[g]
            #  mgals_yz[j, k] += mg[g]
            #  mgals_xy[k, i] += mg[g]


            # CIC interpolation
            #----------------------------------
            #  idown = int((xg[g] - dxhalf)/dx)
            #  jdown = int((yg[g] - dxhalf)/dx)
            #  kdown = int((zg[g] - dxhalf)/dx)
            #  iup = idown + 1
            #  jup = jdown + 1
            #  kup = kdown + 1
            #
            #  #  get volume fractions
            #  xup = xg[g] + dxhalf - iup*dx
            #  yup = yg[g] + dxhalf - jup*dx
            #  zup = zg[g] + dxhalf - kup*dx
            #
            #  if iup == nbins:
            #      iup = idown
            #  if jup == nbins:
            #      jup = jdown
            #  if kup == nbins:
            #      kup = kdown
            #
            #  rho = mg[g] / area
            #
            #  mgals_xy[iup, jup]      += xup      * yup       * rho
            #  mgals_xy[iup, jdown]    += xup      * (dx-yup)  * rho
            #  mgals_xy[idown, jup]    += (dx-xup) * yup       * rho
            #  mgals_xy[idown, jdown]  += (dx-xup) * (dx-yup)  * rho
            #
            #  mgals_yz[jup, kup]      += yup      * zup       * rho
            #  mgals_yz[jup, kdown]    += yup      * (dx-zup)  * rho
            #  mgals_yz[jdown, kup]    += (dx-yup) * zup       * rho
            #  mgals_yz[jdown, kdown]  += (dx-yup) * (dx-zup)  * rho
            #
            #  mgals_zx[kup, iup]      += zup      * xup       * rho
            #  mgals_zx[kup, idown]    += zup      * (dx-xup)  * rho
            #  mgals_zx[kdown, iup]    += (dx-zup) * xup       * rho
            #  mgals_zx[kdown, idown]  += (dx-zup) * (dx-xup)  * rho


            # Fixed galaxy radius
            #----------------------------------

            ilo = int((xg[g] - rgal)//dx)       # round down
            ihi = int((xg[g] + rgal)//dx) + 1   # round up
            dxi = get_volume_fractions(ilo, ihi, xg[g])

            jlo = int((yg[g] - rgal)//dx)       # round down
            jhi = int((yg[g] + rgal)//dx) + 1   # round up
            dyi = get_volume_fractions(jlo, jhi, yg[g])

            klo = int((zg[g] - rgal)//dx)       # round down
            khi = int((zg[g] + rgal)//dx) + 1   # round up
            dzi = get_volume_fractions(klo, khi, zg[g])

            rho = mg[g] / (2*rgal)**2


            totsurf = 0
            for i, xi in dxi:
                for j, yi in dyi:
                    mgals_xy[i, j] += xi * yi * rho
                    totsurf += xi * yi
            for j, yi in dyi:
                for k, zi in dzi:
                    mgals_yz[j, k] += yi * zi * rho
            for k, zi in dzi:
                for i, xi in dxi:
                    mgals_zx[i, j] += zi * xi * rho


        print("Mass check galaxies", mgals_xy.sum()/mg.sum(), mgals_yz.sum()/mg.sum(), mgals_zx.sum()/mg.sum())



        # create orphan mass projections
        morphs_xy = np.zeros((nbins, nbins), dtype='float')
        morphs_yz = np.zeros((nbins, nbins), dtype='float')
        morphs_zx = np.zeros((nbins, nbins), dtype='float')


        for o in range(mo.shape[0]):

            # NGP interpolation
            #--------------------------

            #  i = int(xo[o]/dx)
            #  j = int(yo[o]/dx)
            #  k = int(zo[o]/dx)
            #  morphs_xy[i, j] += mo[o]
            #  morphs_yz[j, k] += mo[o]
            #  morphs_xy[k, i] += mo[o]



            # CIC interpolation
            #--------------------------
            #  idown = int((xo[o] - dxhalf)/dx)
            #  jdown = int((yo[o] - dxhalf)/dx)
            #  kdown = int((zo[o] - dxhalf)/dx)
            #  iup = idown + 1
            #  jup = jdown + 1
            #  kup = kdown + 1
            #
            #  # get volume fractions
            #  xup = xo[o] + dxhalf - iup*dx
            #  yup = yo[o] + dxhalf - jup*dx
            #  zup = zo[o] + dxhalf - kup*dx
            #
            #  if iup == nbins:
            #      iup = idown
            #  if jup == nbins:
            #      jup = jdown
            #  if kup == nbins:
            #      kup = kdown
            #
            #  rho = mo[o] / area
            #
            #  morphs_xy[iup, jup]      += xup      * yup       * rho
            #  morphs_xy[iup, jdown]    += xup      * (dx-yup)  * rho
            #  morphs_xy[idown, jup]    += (dx-xup) * yup       * rho
            #  morphs_xy[idown, jdown]  += (dx-xup) * (dx-yup)  * rho
            #
            #  morphs_yz[jup, kup]      += yup      * zup       * rho
            #  morphs_yz[jup, kdown]    += yup      * (dx-zup)  * rho
            #  morphs_yz[jdown, kup]    += (dx-yup) * zup       * rho
            #  morphs_yz[jdown, kdown]  += (dx-yup) * (dx-zup)  * rho
            #
            #  morphs_zx[kup, iup]      += zup      * xup       * rho
            #  morphs_zx[kup, idown]    += zup      * (dx-xup)  * rho
            #  morphs_zx[kdown, iup]    += (dx-zup) * xup       * rho
            #  morphs_zx[kdown, idown]  += (dx-zup) * (dx-xup)  * rho



            # Fixed galaxy radius
            #----------------------------------

            ilo = int((xo[o] - rgal)//dx)       # round down
            ihi = int((xo[o] + rgal)//dx) + 1   # round up
            dxi = get_volume_fractions(ilo, ihi, xo[o])

            jlo = int((yo[o] - rgal)//dx)       # round down
            jhi = int((yo[o] + rgal)//dx) + 1   # round up
            dyi = get_volume_fractions(jlo, jhi, yo[o])

            klo = int((zo[o] - rgal)//dx)       # round down
            khi = int((zo[o] + rgal)//dx) + 1   # round up
            dzi = get_volume_fractions(klo, khi, zo[o])

            rho = mo[o] / (2*rgal)**2


            totsurf = 0
            for i, xi in dxi:
                for j, yi in dyi:
                    morphs_xy[i, j] += xi * yi * rho
                    totsurf += xi * yi
            for j, yi in dyi:
                for k, zi in dzi:
                    morphs_yz[j, k] += yi * zi * rho
            for k, zi in dzi:
                for i, xi in dxi:
                    morphs_zx[i, j] += zi * xi * rho


        print("Mass check orphans", morphs_xy.sum()/mo.sum(), morphs_yz.sum()/mo.sum(), morphs_zx.sum()/mo.sum())



        # add galaxies to orphans
        morphs_xy += mgals_xy
        morphs_yz += mgals_yz
        morphs_zx += mgals_zx


        # transform projections into surface densities
        area_per_cell = dx**2 / r200**2 # do everything in units of r200

        ngals_xy /= area_per_cell
        ngals_yz /= area_per_cell
        ngals_zx /= area_per_cell
        norphs_xy /= area_per_cell
        norphs_yz /= area_per_cell
        norphs_zx /= area_per_cell

        mgals_xy /= area_per_cell
        mgals_yz /= area_per_cell
        mgals_zx /= area_per_cell
        morphs_xy /= area_per_cell
        morphs_yz /= area_per_cell
        morphs_zx /= area_per_cell



        r = np.zeros(3*nbins**2, dtype=np.float)
        sigma_ng = np.zeros(3*nbins**2, dtype=np.float)
        sigma_no = np.zeros(3*nbins**2, dtype=np.float)
        sigma_mg = np.zeros(3*nbins**2, dtype=np.float)
        sigma_mo = np.zeros(3*nbins**2, dtype=np.float)

        ind = 0
        for i in range(nbins):
            for j in range(nbins):
                xi = xc - (i+0.5)*dx - xlo
                yi = yc - (j+0.5)*dx - ylo
                r[ind] = np.sqrt(xi**2 + yi**2)
                sigma_ng[ind] = ngals_xy[i,j]
                sigma_no[ind] = norphs_xy[i,j]
                sigma_mg[ind] = mgals_xy[i,j]
                sigma_mo[ind] = morphs_yz[i,j]
                ind += 1

        for j in range(nbins):
            for k in range(nbins):
                yi = yc - (j+0.5)*dx - ylo
                zi = zc - (k+0.5)*dx - zlo
                r[ind] = np.sqrt(zi**2 + yi**2)
                sigma_ng[ind] = ngals_yz[j,k]
                sigma_no[ind] = norphs_yz[j,k]
                sigma_mg[ind] = mgals_yz[j,k]
                sigma_mo[ind] = morphs_yz[j,k]
                ind += 1

        for k in range(nbins):
            for i in range(nbins):
                zi = zc - (k+0.5)*dx - zlo
                xi = xc - (i+0.5)*dx - xlo
                r[ind] = np.sqrt(zi**2 + xi**2)
                sigma_ng[ind] = ngals_zx[k,i]
                sigma_no[ind] = norphs_zx[k,i]
                sigma_mg[ind] = mgals_yz[k,i]
                sigma_mo[ind] = morphs_yz[k,i]
                ind += 1


        binarr = np.linspace(1e-2, 3, nbins)
        r /= r200 # do everything in units of r200


        # histogrammize  data now
        nprof_o, edges1 = np.histogram(r, bins=binarr, weights=sigma_no)
        nprof_g, edges2 = np.histogram(r, bins=binarr, weights=sigma_ng)
        mprof_o, edges3 = np.histogram(r, bins=binarr, weights=sigma_mo)
        mprof_g, edges4 = np.histogram(r, bins=binarr, weights=sigma_mg)
        counts, edges = np.histogram(r, bins=binarr, weights=None)
        r_plot = 0.5*(binarr[1:] + binarr[:-1])

        # save histograms
        # The bins have been hardcoded, so every data generated by this script
        # should be on identical bins relative to r200!
        pickle.dump(r_plot, dump)
        pickle.dump(nprof_g, dump)
        pickle.dump(nprof_o, dump)
        pickle.dump(mprof_g, dump)
        pickle.dump(mprof_o, dump)
        pickle.dump(counts, dump)


        mask = counts > 0
        nprof_o[mask] /= counts[mask]
        nprof_g[mask] /= counts[mask]
        mprof_o[mask] /= counts[mask]
        mprof_g[mask] /= counts[mask]




        # get curve fits
        mask = nprof_o > 0
        ydata = nprof_o[mask]
        xdata = r_plot[mask]
        rhoguess = 0.8 * (ydata[0] + ydata[-1])
        opt, cov = curve_fit(sigma_nfw, xdata, ydata, p0 = [2, rhoguess], bounds=(0, [1000, 1e30]), maxfev=100000)
        no_c, no_rho = opt

        mask = mprof_o > 0
        ydata = mprof_o[mask]
        xdata = r_plot[mask]
        rhoguess = 0.8 * (ydata[0] + ydata[-1])
        opt, cov = curve_fit(sigma_nfw, xdata, ydata, p0 = [2, rhoguess], bounds=(0, [1000, 1e30]), maxfev=100000)
        mo_c, mo_rho = opt


        mask = nprof_g > 0
        ydata = nprof_g[mask]
        xdata = r_plot[mask]
        rhoguess = 0.8 * (ydata[0] + ydata[-1])
        opt, cov = curve_fit(sigma_nfw, xdata, ydata, p0 = [2, rhoguess], bounds=(0, [1000, 1e30]), maxfev=100000)
        ng_c, ng_rho = opt

        mask = nprof_g > 0
        ydata = mprof_g[mask]
        xdata = r_plot[mask]
        rhoguess = 0.8 * (ydata[0] + ydata[-1])
        opt, cov = curve_fit(sigma_nfw, xdata, ydata, p0 = [2, rhoguess], bounds=(0, [1000, 1e30]), maxfev=100000)
        mg_c, mg_rho = opt



        # store results
        norphs_by_thresh[I] = nprof_o
        norph_fits[I] = (no_c, no_rho)
        ngals_by_thresh[I] = nprof_g
        ngals_fits[I] = (ng_c, ng_rho)

        morphs_by_thresh[I] = mprof_o
        morph_fits[I] = (mo_c, mo_rho)
        mgals_by_thresh[I] = mprof_g
        mgals_fits[I] = (mg_c, mg_rho)

    dump.close()





    # create plots

    fig = plt.figure(figsize=(8,8))

    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']


    # plot number densities
    ax1 = fig.add_subplot(221)

    for t, thresh in enumerate(threshold_labels):

        nprof_o = norphs_by_thresh[t]
        no_c, no_rho = norph_fits[t]

        mask = nprof_o > 0
        ax1.loglog(r_plot[mask], nprof_o[mask], c=colors[t], alpha=0.4)
        ax1.loglog(r_plot, sigma_nfw(r_plot, no_c, no_rho), "-.", c=colors[t])

    nprof_g = ngals_by_thresh[0]
    ng_c, ng_rho = ngals_fits[0]
    mask = nprof_g > 0
    ax1.loglog(r_plot[mask], nprof_g[mask], c=colors[t+1], alpha=0.4)
    ax1.loglog(r_plot, sigma_nfw(r_plot, ng_c, ng_rho), "-.", c=colors[t+1])

    ax1.set_ylabel(r"$\Sigma_r$ [$r_{200}\ ^{-2}$]")
    ax1.set_title("surface galaxy number density of halo "+halonr)
    ax1.set_xlabel(r'$r/r_{200}$')
    ax1.set_xlim(0.1, 3)



    # make an emptly plot with the legend
    ax2 = fig.add_subplot(222)

    x = [-1, -1]
    y = [-1, -1]
    for t, thresh in enumerate(threshold_labels):
        no_c, no_rho = norph_fits[t]
        ax2.plot(x, y, c=colors[t], alpha=0.4, label="including orphan galaxies, "+thresh)
        ax2.plot(x, y, "-.", c=colors[t], label='including orphans fit, c={0:.3f}, '.format(no_c)+thresh)


    ng_c = ngals_fits[0][0]
    thresh = threshold_labels[0]
    ax2.plot(x, y, c=colors[t+1], alpha=0.4, label="excluding orphan galaxies, "+thresh)
    ax2.plot(x, y, c=colors[t+1], label='excluding orphans fit, c={0:.3f}, '.format(ng_c)+thresh)
    ax2.legend()
    ax2.set_xlim(0,1)
    ax2.set_ylim(0,1)
    ax2.set_axis_off()




    # plot mass densities
    ax3 = fig.add_subplot(223)

    for t, thresh in enumerate(threshold_labels):

        mprof_o = morphs_by_thresh[t]
        mo_c, mo_rho = morph_fits[t]

        mask = mprof_o > 0
        ax3.loglog(r_plot[mask], mprof_o[mask], c=colors[t], alpha=0.4)
        ax3.loglog(r_plot, sigma_nfw(r_plot, mo_c, mo_rho), "-.", c=colors[t])

    mprof_g = mgals_by_thresh[0]
    mg_c, mg_rho = mgals_fits[0]

    mask = mprof_g > 0
    ax3.loglog(r_plot[mask], mprof_g[mask], c=colors[t+1], alpha=0.4)
    ax3.loglog(r_plot, sigma_nfw(r_plot, mg_c, mg_rho), "-.", c=colors[t+1])

    ax3.set_ylabel(r"$\Sigma_r$ [$M_\odot \ r_{200}\ ^{-2}$]")
    ax3.set_title("surface mass density of halo "+halonr)
    ax3.set_xlabel(r'$r/r_{200}$')
    ax3.set_xlim(0.1, 3)




    # make an emptly plot with the legend
    ax4 = fig.add_subplot(224)

    x = [-1, -1]
    y = [-1, -1]
    for t, thresh in enumerate(threshold_labels):
        mo_c, mo_rho = morph_fits[t]
        ax4.plot(x, y, c=colors[t], alpha=0.4, label="including orphan galaxies, "+thresh)
        ax4.plot(x, y, "-.", c=colors[t], label='including orphans fit, c={0:.3f}, '.format(mo_c)+thresh)


    mg_c = mgals_fits[0][0]
    thresh = threshold_labels[0]
    ax4.plot(x, y, c=colors[t+1], alpha=0.4, label="excluding orphan galaxies, "+thresh)
    ax4.plot(x, y, c=colors[t+1], label='excluding orphans fit, c={0:.3f}, '.format(mg_c)+thresh)
    ax4.legend()
    ax4.set_xlim(0,1)
    ax4.set_ylim(0,1)
    ax4.set_axis_off()






    # save image
    fname = "radial-profiles-"+outputnr.zfill(5)+"-halo-"+halonr+".png"
    plt.tight_layout()
    #  plt.savefig(fname)
    #  print("saved figure", fname)
    plt.show()
    plt.close()




    return






def get_profiles_old():

    #  nonlocal r_plot, r_mass
    #  nonlocal mprof_g, mprof_o
    #  nonlocal mgals_av, morphs_av, mgals_av_old, morphs_av_old

    #  for I, mthresh in enumerate(thresholds):
    mthresh = 0

    mask = mg_all > mthresh
    xg = xg_all[mask]
    yg = yg_all[mask]
    zg = zg_all[mask]
    mg = mg_all[mask]

    mask = mo_all > mthresh
    xo = xo_all[mask]
    yo = yo_all[mask]
    zo = zo_all[mask]
    mo = mo_all[mask]


    # create galaxy number projections

    ngals_xy, xedges, yedges = np.histogram2d(xg, yg, bins=nbins, weights=None, range=((xlo, xhi), (ylo, yhi)))
    norphs_xy, xedges, yedges = np.histogram2d(xo, yo, bins=nbins, weights=None, range=((xlo, xhi), (ylo, yhi)))
    norphs_xy += ngals_xy

    ngals_yz, xedges, yedges = np.histogram2d(yg, zg, bins=nbins, weights=None, range=((ylo, yhi), (zlo, zhi)))
    norphs_yz, xedges, yedges = np.histogram2d(yo, zo, bins=nbins, weights=None, range=((ylo, yhi), (zlo, zhi)))
    norphs_yz += ngals_yz

    ngals_zx, xedges, yedges = np.histogram2d(zg, xg, bins=nbins, weights=None, range=((zlo, zhi), (xlo, xhi)))
    norphs_zx, xedges, yedges = np.histogram2d(zo, xo, bins=nbins, weights=None, range=((zlo, zhi), (xlo, xhi)))
    norphs_zx += ngals_zx

    # average projections
    ngals_av = (ngals_xy + ngals_yz + ngals_zx) / 3
    norphs_av = (norphs_xy + norphs_yz + norphs_zx) / 3


    # TODO: temporary
    # create galaxy mass projections
    mgals_xy_old, xedges, yedges = np.histogram2d(xg, yg, bins=nbins, weights=mg, range=((xlo, xhi), (ylo, yhi)))
    morphs_xy_old, xedges, yedges = np.histogram2d(xo, yo, bins=nbins, weights=mo, range=((xlo, xhi), (ylo, yhi)))
    morphs_xy_old += mgals_xy_old

    mgals_yz_old, xedges, yedges = np.histogram2d(yg, zg, bins=nbins, weights=mg, range=((ylo, yhi), (zlo, zhi)))
    morphs_yz_old, xedges, yedges = np.histogram2d(yo, zo, bins=nbins, weights=mo, range=((ylo, yhi), (zlo, zhi)))
    morphs_yz_old += mgals_yz_old

    mgals_zx_old, xedges, yedges = np.histogram2d(zg, xg, bins=nbins, weights=mg, range=((zlo, zhi), (xlo, xhi)))
    morphs_zx_old, xedges, yedges = np.histogram2d(zo, xo, bins=nbins, weights=mo, range=((zlo, zhi), (xlo, xhi)))
    morphs_zx_old += mgals_zx_old

    # average projections
    mgals_av_old = (mgals_xy_old + mgals_yz_old + mgals_zx_old) / 3
    morphs_av_old = (morphs_xy_old + morphs_yz_old + morphs_zx_old) / 3


    #  mtot = mo.sum() + mg.sum()
    #  mtot = mo.sum()
    #  print("old check", morphs_xy_old.sum()/mtot, morphs_yz_old.sum()/mtot, morphs_zx_old.sum()/mtot, morphs_av_old.sum()/mtot)




    # transform projections into surface densities
    area_per_cell = dx**2 / r200**2 # do everything in units of r200
    ngals_av /= area_per_cell
    norphs_av /= area_per_cell
    mgals_av /= area_per_cell
    morphs_av /= area_per_cell
    mgals_av_old /= area_per_cell
    morphs_av_old /= area_per_cell


    # get all surface densities in a 1D array along with a distance array

    r = np.zeros(nbins**2, dtype=np.float)
    sigma_ng = np.zeros(nbins**2, dtype=np.float)
    sigma_no = np.zeros(nbins**2, dtype=np.float)
    sigma_mg = np.zeros(nbins**2, dtype=np.float)
    sigma_mo = np.zeros(nbins**2, dtype=np.float)

    ind = 0
    for i in range(nbins):
        for j in range(nbins):
            xi = xc - ((i+0.5)*dx + xlo)
            yi = yc - ((j+0.5)*dx + ylo)
            r[ind] = np.sqrt(xi**2 + yi**2)
            sigma_ng[ind] = ngals_av[i,j]
            sigma_no[ind] = norphs_av[i,j]
            sigma_mg[ind] = mgals_av_old[i,j]
            sigma_mo[ind] = morphs_av_old[i,j]
            ind += 1

    print("sum sigma", sigma_mo.sum()/mo.sum())

    binarr = np.linspace(0., 4, nbins)
    #  binarr = np.linspace(0., rmax_cells-0.5*dx, nbins+1) / r200
    r /= r200 # do everything in units of r200


    # histogrammize  data now
    nprof_o, edges1 = np.histogram(r, bins=binarr, weights=sigma_no)
    nprof_g, edges2 = np.histogram(r, bins=binarr, weights=sigma_ng)
    mprof_o, edges3 = np.histogram(r, bins=binarr, weights=sigma_mo)
    mprof_g, edges4 = np.histogram(r, bins=binarr, weights=sigma_mg)
    counts, edges = np.histogram(r, bins=binarr, weights=None)
    r_plot = 0.5*(binarr[1:] + binarr[:-1])

    print("Total mass in mprof")

    # save histograms
    # The bins have been hardcoded, so every data generated by this script
    # should be on identical bins relative to r200!
    #  pickle.dump(r_plot, dump)
    #  pickle.dump(nprof_g, dump)
    #  pickle.dump(nprof_o, dump)
    #  pickle.dump(mprof_g, dump)
    #  pickle.dump(mprof_o, dump)
    #  pickle.dump(counts, dump)



    # normalize histograms such that you get average values, not summed up ones
    mask = counts > 0
    nprof_o = nprof_o[mask] / counts[mask]
    nprof_g = nprof_g[mask] / counts[mask]
    mprof_o = mprof_o[mask] / counts[mask]
    mprof_g = mprof_g[mask] / counts[mask]
    r_plot = r_plot[mask]

    return

