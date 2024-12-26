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
# plot_radial_profiles_with_std.py 0067 123456
#===========================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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

print_table_data = True
table_data_r200_cutoff = True
make_plots = False # skip actual plotting if you only want output for stacked profiles
# use also upper mass thresholds for galaxies
use_upper_threshold = True


nlogbins = 40
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
def fit_nfw(xdata, ydata, std):
#================================
    """
    Get a fit to NFW surface density profile
    xdata: x values
    ydata: y values
    std: standard deviation of y values

    Returns: 
        (c, rho): Parameters to be fitted
    """
    
    cmax = 20.
    rhoguess = 0.5 * (ydata[0] + ydata[-1])
    opt, cov = curve_fit(
                            sigma_nfw, 
                            xdata, 
                            ydata, 
                            sigma=std, 
                            p0 = [2., rhoguess], 
                            bounds=(0, [cmax, 1e30]), 
                            maxfev=10000000
                        )
    c, rho = opt
    if abs(c/cmax - 1) < 0.01:
        # try again with reduced dataset
        if xdata.shape[0] > 6:
            print("restarting fit")
            if std is None:
                return fit_nfw(xdata[1:-1], ydata[1:-1], std)
            else:
                return fit_nfw(xdata[1:-1], ydata[1:-1], std[1:-1])
    return c, rho


#================================================
def get_hists(rxy, ryz, rzx, weights, bins):
#================================================

    xy, bin_edges = np.histogram(rxy, bins=bins, weights=weights) 
    yz, bin_edges = np.histogram(ryz, bins=bins, weights=weights) 
    zx, bin_edges = np.histogram(rzx, bins=bins, weights=weights) 

    surfaces = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)
    bin_centres = (bin_edges[:-1] + bin_edges[1:]) * 0.5

    xy = xy / surfaces
    yz = yz / surfaces
    zx = zx / surfaces

    return xy, yz, zx, bin_centres



#================================================
def get_mean_and_errors(xydata, yzdata, zxdata):
#================================================

    stack = np.vstack((xydata, yzdata, zxdata))
    mean = np.mean(stack, axis=0)
    std = np.std(stack, axis=0)
    return mean, std

#======================
def getfloat(f):
#======================
    line = f.readline()
    line = line.strip()
    return float(line)


#======================================================================
def plot_radial_profile(outputnr, halonr):
#======================================================================
    """
    Plots radial profile.
    Create a simple radial histogram.
    """

    nr = str(outputnr).zfill(5)


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
    rg_xy = (xg - xc)**2 + (yg - yc)**2
    rg_yz = (yg - yc)**2 + (zg - zc)**2
    rg_zx = (zg - zc)**2 + (xg - xc)**2
    ro_xy = (xo - xc)**2 + (yo - yc)**2
    ro_yz = (yo - yc)**2 + (zo - zc)**2
    ro_zx = (zo - zc)**2 + (xo - xc)**2

    if table_data_r200_cutoff:
        rg = np.sqrt((xg - xc)**2 + (yg - yc)**2 + (zg - zc)**2)/r200
        ro = np.sqrt((xo - xc)**2 + (yo - yc)**2 + (zo - zc)**2)/r200

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
    rg_xy = np.sqrt(rg_xy) / r200
    rg_yz = np.sqrt(rg_yz) / r200
    rg_zx = np.sqrt(rg_zx) / r200
    ro_xy = np.sqrt(ro_xy) / r200
    ro_yz = np.sqrt(ro_yz) / r200
    ro_zx = np.sqrt(ro_zx) / r200

    # set limits by hand
    rmin = 0.02
    rmax = 1.

    #  print("rmin in kpc:", rmin * r200 * 1000)

    #  print("r_min:", rmin, "r200")
    #  print("r_min:", rmin*r200*1000, "kpc") # default unit is Mpc

    logbins = np.logspace(np.log10(rmin), np.log10(rmax*2), num=nlogbins)
    # use varying bin widths to give more emphasis on outer regions of a halo
    # number density seems fine. Look to possibly improve mass density!
    #  linbins_number_density = np.hstack((np.linspace(rmin*0.99, 0.1, 10), np.linspace(0.11, rmax*1.01, 40)))
    #  linbins_mass_density = np.hstack((np.linspace(rmin*0.99, 0.2, 5), np.linspace(0.21, rmax*1.01, 100)))
    #  linbins_number_density = np.hstack((np.linspace(rmin*0.99, 0.2, 20), np.linspace(0.21, 0.5, 20), np.linspace(0.51, rmax*1.01, 40)))
    #  linbins_number_density = np.hstack((np.linspace(rmin*0.99, 0.2, 5), np.linspace(0.21, 0.91, 20), np.linspace(0.92, rmax*1.01, 20)))
    linbins_number_density = np.hstack((np.linspace(rmin*0.99, 0.1, 10), np.linspace(0.11, 0.5, 30), np.linspace(0.51, rmax*1.01, 20)))
    linbins_mass_density = np.hstack((np.linspace(rmin*0.99, 0.1, 10), np.linspace(0.11, 0.5, 30), np.linspace(0.51, rmax*1.01, 20)))

    surface_densities_gxy = []
    surface_densities_gyz = []
    surface_densities_gzx = []
    surface_densities_gmean = []
    surface_densities_gstd = []

    surface_densities_oxy = []
    surface_densities_oyz = []
    surface_densities_ozx = []
    surface_densities_omean = []
    surface_densities_ostd = []
    
    number_densities_gxy = []
    number_densities_gyz = []
    number_densities_gzx = []
    number_densities_gmean = []
    number_densities_gstd = []
   
    number_densities_oxy = []
    number_densities_oyz = []
    number_densities_ozx = []
    number_densities_omean = []
    number_densities_ostd = []

    lin_surface_densities_gxy = []
    lin_surface_densities_gyz = []
    lin_surface_densities_gzx = []
    lin_surface_densities_gmean = []
    lin_surface_densities_gstd = []

    lin_surface_densities_oxy = []
    lin_surface_densities_oyz = []
    lin_surface_densities_ozx = []
    lin_surface_densities_omean = []
    lin_surface_densities_ostd = []

    lin_number_densities_gxy = []
    lin_number_densities_gyz = []
    lin_number_densities_gzx = []
    lin_number_densities_gmean = []
    lin_number_densities_gstd = []

    lin_number_densities_oxy = []
    lin_number_densities_oyz = []
    lin_number_densities_ozx = []
    lin_number_densities_omean = []
    lin_number_densities_ostd = []


    surface_density_fit_g = []
    surface_density_fit_o = []
    number_density_fit_g = []
    number_density_fit_o = []


    thresholds = [0, 1e9, 1e10]
    if use_upper_threshold:
        threshold_labels = [r"$M_* < 10^9 M_\odot$", r"$10^9 M_\odot < M_* < 10^{10} M_\odot$", r"$M_* > 10^{10} M_\odot$"]
    else:
        threshold_labels = ["no threshold", r"$M_* > 10^9 M_\odot$", r"$M_* > 10^{10} M_\odot$"]


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






    for tind, t in enumerate(thresholds):
        #------------------------
        # Histogram Data
        #------------------------

        if use_upper_threshold:
            if tind == 0:
                choice_g = mg < thresholds[tind+1] # we need to have at least 2 thresholds!
                choice_o = mo < thresholds[tind+1]
            elif tind == len(thresholds) - 1:
                choice_g = mg > t
                choice_o = mo > t
            else:
                choice_g = np.logical_and((mg > t), (mg < thresholds[tind+1]))
                choice_o = np.logical_and((mo > t), (mo < thresholds[tind+1]))
        else:
            choice_g = mg > t
            choice_o = mo > t
            
        #  print("Threshold", t, "keeping",
        #          np.count_nonzero(choice_g), "/", choice_g.shape[0], "galaxies;",
        #          np.count_nonzero(choice_o), "/", choice_o.shape[0], "orphans;")

        number_gxy, number_gyz, number_gzx, bin_centres = get_hists(rg_xy[choice_g], rg_yz[choice_g], rg_zx[choice_g], None, logbins)
        mass_gxy, mass_gyz, mass_gzx, bin_edges = get_hists(rg_xy[choice_g], rg_yz[choice_g], rg_zx[choice_g], mg[choice_g], logbins)

        number_oxy, number_oyz, number_ozx, bin_centres = get_hists(ro_xy[choice_o], ro_yz[choice_o], ro_zx[choice_o], None, logbins)
        mass_oxy, mass_oyz, mass_ozx, bin_edges = get_hists(ro_xy[choice_o], ro_yz[choice_o], ro_zx[choice_o], mo[choice_o], logbins)

        # add galaxies to orphans
        mass_oxy += mass_gxy
        mass_oyz += mass_gyz
        mass_ozx += mass_gzx
        number_oxy += number_gxy
        number_oyz += number_gyz
        number_ozx += number_gzx

        number_g_mean, number_g_std = get_mean_and_errors(number_gxy, number_gyz, number_gzx)
        mass_g_mean, mass_g_std = get_mean_and_errors(mass_gxy, mass_gyz, mass_gzx)

        number_o_mean, number_o_std = get_mean_and_errors(number_oxy, number_oyz, number_ozx)
        mass_o_mean, mass_o_std = get_mean_and_errors(mass_oxy, mass_oyz, mass_ozx)

        surface_densities_gxy.append(mass_gxy)
        surface_densities_gyz.append(mass_gyz)
        surface_densities_gzx.append(mass_gzx)
        surface_densities_gmean.append(mass_g_mean)
        surface_densities_gstd.append(mass_g_std)

        surface_densities_oxy.append(mass_oxy)
        surface_densities_oyz.append(mass_oyz)
        surface_densities_ozx.append(mass_ozx)
        surface_densities_omean.append(mass_o_mean)
        surface_densities_ostd.append(mass_o_std)
        
        number_densities_gxy.append(number_gxy)
        number_densities_gyz.append(number_gyz)
        number_densities_gzx.append(number_gzx)
        number_densities_gmean.append(number_g_mean)
        number_densities_gstd.append(number_g_std)
       
        number_densities_oxy.append(number_oxy)
        number_densities_oyz.append(number_oyz)
        number_densities_ozx.append(number_ozx)
        number_densities_omean.append(number_o_mean)
        number_densities_ostd.append(number_o_std)



        #------------------------------------
        # Histogram data lineraly as well
        #------------------------------------

        lin_number_gxy, lin_number_gyz, lin_number_gzx, lin_bin_number_centres = get_hists(rg_xy[choice_g], rg_yz[choice_g], rg_zx[choice_g], None, linbins_number_density)
        lin_mass_gxy, lin_mass_gyz, lin_mass_gzx, lin_bin_mass_centres = get_hists(rg_xy[choice_g], rg_yz[choice_g], rg_zx[choice_g], mg[choice_g], linbins_mass_density)

        lin_number_oxy, lin_number_oyz, lin_number_ozx, lin_bin_number_centres = get_hists(ro_xy[choice_o], ro_yz[choice_o], ro_zx[choice_o], None, linbins_number_density)
        lin_mass_oxy, lin_mass_oyz, lin_mass_ozx, lin_bin_mass_centres = get_hists(ro_xy[choice_o], ro_yz[choice_o], ro_zx[choice_o], mo[choice_o], linbins_mass_density)

        # add galaxies to orphans
        lin_mass_oxy += lin_mass_gxy
        lin_mass_oyz += lin_mass_gyz
        lin_mass_ozx += lin_mass_gzx
        lin_number_oxy += lin_number_gxy
        lin_number_oyz += lin_number_gyz
        lin_number_ozx += lin_number_gzx

        lin_number_g_mean, lin_number_g_std = get_mean_and_errors(lin_number_gxy, lin_number_gyz, lin_number_gzx)
        lin_mass_g_mean, lin_mass_g_std = get_mean_and_errors(lin_mass_gxy, lin_mass_gyz, lin_mass_gzx)

        lin_number_o_mean, lin_number_o_std = get_mean_and_errors(lin_number_oxy, lin_number_oyz, lin_number_ozx)
        lin_mass_o_mean, lin_mass_o_std = get_mean_and_errors(lin_mass_oxy, lin_mass_oyz, lin_mass_ozx)

        lin_surface_densities_gmean.append(lin_mass_g_mean)
        lin_surface_densities_gstd.append(lin_mass_g_std)

        lin_surface_densities_omean.append(lin_mass_o_mean)
        lin_surface_densities_ostd.append(lin_mass_o_std)
        
        lin_number_densities_gmean.append(lin_number_g_mean)
        lin_number_densities_gstd.append(lin_number_g_std)
       
        lin_number_densities_omean.append(lin_number_o_mean)
        lin_number_densities_ostd.append(lin_number_o_std)


        lin_surface_densities_gxy.append(lin_mass_gxy)
        lin_surface_densities_gyz.append(lin_mass_gyz)
        lin_surface_densities_gzx.append(lin_mass_gzx)

        lin_surface_densities_oxy.append(lin_mass_oxy)
        lin_surface_densities_oyz.append(lin_mass_oyz)
        lin_surface_densities_ozx.append(lin_mass_ozx)
        
        lin_number_densities_gxy.append(lin_number_gxy)
        lin_number_densities_gyz.append(lin_number_gyz)
        lin_number_densities_gzx.append(lin_number_gzx)
       
        lin_number_densities_oxy.append(lin_number_oxy)
        lin_number_densities_oyz.append(lin_number_oyz)
        lin_number_densities_ozx.append(lin_number_ozx)




        #------------------------
        # Get Fits
        #------------------------
        
        if make_plots:
            mask = lin_number_o_mean > 0
            ydata = lin_number_o_mean[mask]
            xdata = lin_bin_number_centres[mask]
            #  std = lin_number_o_std[mask]
            #  stdmin = std[std > 0.].min()
            #  std[std < stdmin] = 1e-3*stdmin
            std = None
            number_density_fit_o.append(fit_nfw(xdata, ydata, std))

            mask = lin_mass_o_mean > 0
            ydata = lin_mass_o_mean[mask]
            xdata = lin_bin_mass_centres[mask]
            #  std = lin_mass_o_std[mask]
            #  stdmin = std[std > 0.].min()
            #  std[std < stdmin] = 1e-3*stdmin
            std = None
            surface_density_fit_o.append(fit_nfw(xdata, ydata, std))

            mask = lin_number_g_mean > 0
            ydata = lin_number_g_mean[mask]
            xdata = lin_bin_number_centres[mask]
            #  std = lin_number_g_std[mask]
            #  stdmin = std[std > 0.].min()
            #  std[std < stdmin] = 1e-3*stdmin
            std = None
            number_density_fit_g.append(fit_nfw(xdata, ydata, std))

            mask = lin_mass_g_mean > 0
            ydata = lin_mass_g_mean[mask]
            xdata = lin_bin_mass_centres[mask]
            #  std = lin_mass_g_std[mask]
            #  stdmin = std[std > 0.].min()
            #  std[std < stdmin] = 1e-3*stdmin
            std = None
            surface_density_fit_g.append(fit_nfw(xdata, ydata, std))



    # open dump file
    fname = "radial-profiles-data-for-stacked-profile-std-"+outputnr.zfill(5)+"-halo-"+halonr+".pkl"
    dump = open(fname, 'wb')

    pickle.dump(thresholds, dump)
    pickle.dump(threshold_labels, dump)

    pickle.dump(surface_densities_gxy, dump)
    pickle.dump(surface_densities_gyz, dump)
    pickle.dump(surface_densities_gzx, dump)

    pickle.dump(surface_densities_oxy, dump)
    pickle.dump(surface_densities_oyz, dump)
    pickle.dump(surface_densities_ozx, dump)

    pickle.dump(number_densities_gxy, dump)
    pickle.dump(number_densities_gyz, dump)
    pickle.dump(number_densities_gzx, dump)
   
    pickle.dump(number_densities_oxy, dump)
    pickle.dump(number_densities_oyz, dump)
    pickle.dump(number_densities_ozx, dump)

    pickle.dump(lin_surface_densities_gxy, dump)
    pickle.dump(lin_surface_densities_gyz, dump)
    pickle.dump(lin_surface_densities_gzx, dump)

    pickle.dump(lin_surface_densities_oxy, dump)
    pickle.dump(lin_surface_densities_oyz, dump)
    pickle.dump(lin_surface_densities_ozx, dump)

    pickle.dump(lin_number_densities_gxy, dump)
    pickle.dump(lin_number_densities_gyz, dump)
    pickle.dump(lin_number_densities_gzx, dump)

    pickle.dump(lin_number_densities_oxy, dump)
    pickle.dump(lin_number_densities_oyz, dump)
    pickle.dump(lin_number_densities_ozx, dump)

    pickle.dump(bin_centres, dump)
    pickle.dump(lin_bin_number_centres, dump)
    pickle.dump(lin_bin_mass_centres, dump)

    pickle.dump(use_upper_threshold, dump)

    dump.close()








    #----------------------
    # create plots
    #----------------------

    if make_plots:

        fig = plt.figure(figsize=(10,10))

        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']


        # plot number densities
        ax1 = fig.add_subplot(221)
        ax1.set_title("Surface Number Densities of Galaxies with Hosts")
        ax1.set_ylabel(r"$\Sigma_n$ [$r_{200}\ ^{-2}$]")

        for i, t in enumerate(thresholds):
            number_density_g = number_densities_gmean[i]
            err = number_densities_gstd[i]
            mask = number_density_g >= 0
            ax1.errorbar(
                bin_centres[mask], 
                number_density_g[mask], 
                yerr=err[mask], 
                label="galaxies with hosts "+threshold_labels[i], 
                c = colors[i], 
                alpha = alpha, 
                )

            for data in number_densities_gxy, number_densities_gyz, number_densities_gzx:
                mask = data[i] > 0
                ax1.loglog(
                            bin_centres[mask],
                            data[i][mask],
                            c=colors[i],
                            alpha = 0.2 * alpha
                        )

            c, rho = number_density_fit_g[i]
            ax1.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                lw=3, 
                )

        # plot number densities
        ax2 = fig.add_subplot(222)
        ax2.set_title("Surface Number Densities of Galaxies with Hosts and Orphans")
        ax2.set_ylabel(r"$\Sigma_n$ [$r_{200}\ ^{-2}$]")

        for i, t in enumerate(thresholds):
            number_density_o = number_densities_omean[i]
            err = number_densities_ostd[i]
            mask = number_density_o >= 0
            ax2.errorbar(
                        bin_centres[mask], 
                        number_density_o[mask], 
                        yerr = err[mask], 
                        label="all galaxies, including orphans "+threshold_labels[i],
                        c = colors[i],
                        alpha = alpha, 
                        )

            for data in number_densities_oxy, number_densities_oyz, number_densities_ozx:
                mask = data[i] > 0
                ax2.loglog(
                            bin_centres[mask],
                            data[i][mask],
                            c=colors[i],
                            alpha = 0.2 * alpha
                        )


            c, rho = number_density_fit_o[i]
            ax2.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                lw=3, 
                )


        # make an emptly plot with the legend
        ax3 = fig.add_subplot(223)
        ax3.set_title("Surface Mass Densities of Galaxies with Hosts")
        ax3.set_ylabel(r"$\Sigma_m$ [$M_\odot$ $r_{200}\ ^{-2}$]")
        for i, t in enumerate(thresholds):
            surface_density_g = surface_densities_gmean[i]
            err = surface_densities_gstd[i]
            mask = surface_density_g > 0
            ax3.errorbar(
                        bin_centres[mask], 
                        surface_density_g[mask], 
                        yerr=err[mask], 
                        label="galaxies with hosts "+threshold_labels[i], 
                        c = colors[i], 
                        alpha = alpha, 
                        )

            for data in surface_densities_gxy, surface_densities_gyz, surface_densities_gzx:
                mask = data[i] > 0
                ax3.loglog(
                            bin_centres[mask],
                            data[i][mask],
                            c=colors[i],
                            alpha = 0.2 * alpha
                        )

            c, rho = surface_density_fit_g[i]
            ax3.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                lw=3, 
                )


        ax4 = fig.add_subplot(224)
        ax4.set_title("Surface Mass Densities of Galaxies with Hosts and orphans")
        ax4.set_ylabel(r"$\Sigma_m$ [$M_\odot$ $r_{200}\ ^{-2}$]")
        for i, t in enumerate(thresholds):
            surface_density_o = surface_densities_omean[i]
            err = surface_densities_ostd[i]
            mask = surface_density_o > 0

            ax4.errorbar(
                bin_centres[mask], 
                surface_density_o[mask], 
                yerr=err[mask], 
                label="including orphans "+threshold_labels[i], 
                c = colors[i], 
                alpha = alpha, 
                )

            for data in surface_densities_oxy, surface_densities_oyz, surface_densities_ozx:
                mask = data[i] > 0
                ax4.loglog(
                            bin_centres[mask],
                            data[i][mask],
                            c=colors[i],
                            alpha = 0.2 * alpha
                        )

            c, rho = surface_density_fit_o[i]
            ax4.loglog(
                bin_centres, 
                sigma_nfw(bin_centres, c, rho), 
                "-.",
                c=colors[i],
                label="fit with c={0:.2f}".format(c),
                alpha = fitalpha, 
                lw=3, 
                )


        for ax in fig.axes:
            ax.legend()
            ax.set_xlabel(r"$r/r_{200}$")
            ax.grid()




        # save image
        fname = "radial-profiles-std-"+outputnr.zfill(5)+"-halo-"+halonr+".png"
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
