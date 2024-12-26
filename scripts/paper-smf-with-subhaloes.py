#!/usr/bin/env python3
#---------------------------------------------------
# Plots results obtained from get_smf.f03 for the
# report. This script looks for the following
# directories:
#   galaxy-512-new-cosmo-69MPC
#   galaxy-512-new-cosmo-100MPC
#---------------------------------------------------

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()
import os

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                                'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)

# for legend
from matplotlib.font_manager import FontProperties # for legend
fontP=FontProperties()
fontP.set_size('x-small')
fontP.set_family('serif')

Mpc = 3.086e24 # Mpc/cm
h_cosmo = 0.673

plot_obs_only = False
#  plot_obs_only = True




redshift_chosen = 0.1 #, 0.25, 0.45, 0.7, 0.9, 1.25, 1.65, 2.25 ] # list of redshift bin centers
snapshot = 61 #,  53, 44, 36, 31, 24, 19, 14]   # snapshot numbers corresponding to chosen redshifts
simdirs = ['galaxy-512-new-cosmo-69MPC', 'galaxy-512-new-cosmo-100MPC'] # root directories of simulations
simname = {} # dict for axis labelling
simname[simdirs[0]] = 'G69'
simname[simdirs[1]] = 'G100'




obsdir = '/home/mivkov/UZH/Masterarbeit/masterarbeit/files/observational_data'

allgalsobsfile = 'moustakasSMF_z0.1.txt'
allgalsobslabel = 'Moustakas+13'
#  allgalsobslabel = 'https://arxiv.org/pdf/1301.1688.pdf Table 3'

centralsobsfile = 'behroozi-2013-data-compilation/smf_ms/moustakas_z0.105.smf-logarithmic'
centralsobslabel = "Behroozi+13 data for Moustakas+13"

yangobsfile = 'yang08_SMFs.txt'
yangobslabel = 'Yang+08'



#==================================
def get_my_smf(simdir):
#==================================
    """
    Read in my SMF from smf.txt
    """

    snapdir = os.path.join(simdir, 'output_'+str(snapshot).zfill(5))
    smffile = os.path.join(snapdir, 'smf-including-satellites.txt')
    #  smffile = os.path.join(snapdir, 'smf-new.txt')
    infofname = os.path.join(snapdir, 'info_'+str(snapshot).zfill(5)+'.txt')

    #------------------------------
    # first read infofile
    #------------------------------

    infofile = open(infofname)
    for j in range(9):
        infofile.readline() # skip first 9 lines

    # get expansion factor
    aline = infofile.readline()
    astring, equal, aval = aline.partition("=")
    afloat = float(aval)
    a_exp = afloat

    for j in range(5):
        infofile.readline() # skip first 9 lines

    lline = infofile.readline()
    lstring, equal, lval=lline.partition('=')
    lfloat = float(lval)
    unit_l = lfloat/Mpc

    infofile.close()

    redshift = 1.0/a_exp - 1.0


    #-----------------------------
    # Now read SMF
    #-----------------------------

    #  masses, smf_central, smf_satellite, smf_orph = np.loadtxt(smffile, dtype='float', skiprows=1, usecols=([0, 1, 2, 3]), unpack=True)
    masses, smf_central, smf_satellite, smf_orph, smf_orph_filtered = np.loadtxt(smffile, dtype='float', skiprows=1, usecols=([0, 1, 2, 3, 4]), unpack=True)

    print("Read in from:", smffile)
    print("   Central:", int(smf_central.sum()))
    print("   Satellite: ", int(smf_satellite.sum()))
    print("   Haloes:", int(smf_central.sum() + smf_satellite.sum()))
    print("   Orphans:", int(smf_orph.sum()))
    print("   Filtered Orphans:", int(smf_orph_filtered.sum()))
    print("   Tot:", int(smf_central.sum() + smf_satellite.sum() + smf_orph.sum() ))


    #-----------------------------
    # Process data
    #-----------------------------

    masses_plot = 0.5*(masses[1:]+masses[:-1])
    dex_divide = masses[1:]-masses[:-1]

    # Note: smf are conventionally defined over comoving volume
    vol = (unit_l*(1+redshift))**3

    # for Moustakas 13, don't use 1+z
    #  vol = (unit_l)**3

    smf_central = smf_central[1:].copy() # skip first empty bin
    mask = smf_central > 0
    dd = dex_divide[mask]
    phi_centrals = np.log10(smf_central[mask]/vol/dd)
    masses_centrals = masses_plot[mask]

    smf_satellite = smf_satellite[1:].copy()
    mask = smf_satellite > 0
    dd = dex_divide[mask]
    phi_satellite = np.log10(smf_satellite[mask]/vol/dd)
    masses_satellite = masses_plot[mask]

    smf_all = smf_central.copy() + smf_satellite.copy()
    mask = smf_all > 0
    dd = dex_divide[mask]
    phi_all = np.log10(smf_all[mask]/vol/dd)
    masses_all = masses_plot[mask]


    smf_all_orph = smf_orph[1:].copy()
    smf_all_orph += smf_all
    mask = smf_all_orph > 0
    dd = dex_divide[mask]
    phi_all_orph = np.log10(smf_all_orph[mask]/vol/dd)
    masses_all_orph = masses_plot[mask]

    smf_sat_orph = smf_orph[1:].copy()
    smf_sat_orph += smf_all
    mask = smf_sat_orph > 0
    dd = dex_divide[mask]
    phi_sat_orph = np.log10(smf_sat_orph[mask]/vol/dd)
    masses_sat_orph = masses_plot[mask]
    

    smf_all_orph_filtered = smf_orph_filtered[1:].copy()
    smf_all_orph_filtered += smf_all
    mask = smf_all_orph_filtered > 0
    dd = dex_divide[mask]
    phi_all_orph_filtered = np.log10(smf_all_orph_filtered[mask]/vol/dd)
    masses_all_orph_filtered = masses_plot[mask]

    smf_sat_orph_filtered = smf_orph_filtered[1:].copy()
    smf_sat_orph_filtered += smf_satellite
    mask = smf_sat_orph_filtered > 0
    dd = dex_divide[mask]
    phi_sat_orph_filtered = np.log10(smf_sat_orph_filtered[mask]/vol/dd)
    masses_sat_orph_filtered = masses_plot[mask]
    
    
    return  (masses_all, phi_all), \
            (masses_centrals, phi_centrals),\
            (masses_satellite, phi_satellite), \
            (masses_all_orph, phi_all_orph), \
            (masses_sat_orph, phi_sat_orph), \
            (masses_all_orph_filtered, phi_all_orph_filtered), \
            (masses_sat_orph_filtered, phi_sat_orph_filtered)




#==============
def main():
#==============

    nrows = 1
    ncols = 3

    fig = plt.figure(figsize=(15, 5))
    #  fig = plt.figure(figsize=(16, 8))
    for i in range(nrows*ncols):
        fig.add_subplot(nrows, ncols, i+1)

    ax1, ax2, ax3 = fig.axes


    # plot observations: Yang08 everything
    #----------------------------------------

    obsfile = os.path.join(obsdir, yangobsfile)
    mass, phi_all, phi_all_error, phi_centrals, phi_centrals_error, phi_satellites, phi_satellites_error = np.loadtxt(obsfile, comments='#', unpack=True, dtype=np.float) 

    # for Yang et al 2013
    if plot_obs_only:
        phi_fact = 1e-2
        logphi_all = np.log10(phi_all * phi_fact)
        # computing errors for log plots:
        # assume you're plotting new quantity z = log10(y(x)), y(x) function of x which you use on the x axis
        # then dz \approx d[log10(y)] = d[ ln(y) / ln(10) ] = 1/ln(10) d[ln(y)] = 1/2.303 * dy/y = 0.434 * dy/y
        logphi_all_error = 0.434 * phi_all_error / phi_all
        logphi_centrals = np.log10(phi_centrals * phi_fact)
        logphi_centrals_error = 0.434 * phi_centrals_error / phi_centrals
        logphi_satellites = np.log10(phi_satellites * phi_fact)
        logphi_satellites_error = 0.434 * phi_satellites_error / phi_satellites

    else:
        phi_fact = 1e-2 * h_cosmo**3
        logphi_all = np.log10(phi_all * phi_fact)
        # computing errors for log plots:
        # assume you're plotting new quantity z = log10(y(x)), y(x) function of x which you use on the x axis
        # then dz \approx d[log10(y)] = d[ ln(y) / ln(10) ] = 1/ln(10) d[ln(y)] = 1/2.303 * dy/y = 0.434 * dy/y
        logphi_all_error = 0.434 * phi_all_error / phi_all
        logphi_centrals = np.log10(phi_centrals * phi_fact)
        logphi_centrals_error = 0.434 * phi_centrals_error / phi_centrals
        logphi_satellites = np.log10(phi_satellites * phi_fact)
        logphi_satellites_error = 0.434 * phi_satellites_error / phi_satellites
        mass -= np.log10(h_cosmo**2)

    ax1.errorbar(mass, logphi_all, yerr = logphi_all_error, label=yangobslabel+" all", zorder=0, capsize=2)
    ax2.errorbar(mass, logphi_centrals, yerr = logphi_centrals_error, label=yangobslabel+" centrals", zorder=0, capsize=2)
    ax3.errorbar(mass, logphi_satellites, yerr = logphi_satellites_error, label=yangobslabel+" satellites", zorder=0, capsize=2)




    # plot observations: Moustakas Centrals used in Behroozi+13
    #--------------------------------------------------------------

    obsfile = os.path.join(obsdir, centralsobsfile)
    mass, phi, yerrdown, yerrup = np.loadtxt(obsfile, comments='#', unpack=True, dtype=np.float)

    if plot_obs_only:
        # add h dependency so the units fit the other two observations
        #  M  = M h^2 h^-2 = [M h^2] h^-2
        mass += 2*np.log10(0.7)
        # phi = (phi h^-3) h^3
        phi -= 3*np.log10(0.7)

    #  TODO: replace to only where needed
    ax1.errorbar(mass, phi, yerr = (yerrup, yerrdown), label=centralsobslabel+" centrals;\n plotted for comparison", zorder=0, capsize=2)
    ax2.errorbar(mass, phi, yerr = (yerrup, yerrdown), label=centralsobslabel, zorder=0, capsize=2)



    # plot observations: Moustakas All galaxies 
    #--------------------------------------------

    obsfile = os.path.join(obsdir, allgalsobsfile)
    mass, phi, yerrdown, yerrup = np.loadtxt(obsfile, comments='#', unpack=True, dtype=np.float) 

    if plot_obs_only:
        # mass = (phi h^-3) h^3
        mass_plot = mass + 2*np.log10(0.7)
        # phi = (phi h^-3) h^3
        phi_plot = phi - 3*np.log10(0.7)
        
        # computing errors for log plots:
        # assume you're plotting new quantity z = log10(y(x)), y(x) function of x which you use on the x axis
        # then dz \approx d[log10(y)] = d[ ln(y) / ln(10) ] = 1/ln(10) d[ln(y)] = 1/2.303 * dy/y = 0.434 * dy/y
        # => dy = ln(10) * y * dz; In this case, I only have dz given in the yerrdown, yerrup arrays
        dyup = 2.303 * 10**phi * yerrup
        zerrorup = 0.434 * dyup / 10**phi_plot
        dydown = 2.303 * 10**phi * yerrdown
        zerrordown = 0.434 * dydown / 10**phi_plot

    else:
        # M = (M h^-2) h^2 = M' h^-2; To get M', divide by h^-2
        # mass array contains M(h=0.70)
        # to convert to my h, we first need M h_70^-2
        mass_hminus2 = mass + 2*np.log10(0.7)
        mass_plot = mass_hminus2 - 2*np.log10(h_cosmo)
        # phi = (phi h^-3) h^3
        phi_hminus3 = phi - 3*np.log10(0.7)
        phi_plot = phi_hminus3 + 3*np.log10(h_cosmo)
        
        # computing errors for log plots:
        # assume you're plotting new quantity z = log10(y(x)), y(x) function of x which you use on the x axis
        # then dz \approx d[log10(y)] = d[ ln(y) / ln(10) ] = 1/ln(10) d[ln(y)] = 1/2.303 * dy/y = 0.434 * dy/y
        # => dy = ln(10) * y * dz; In this case, I only have dz given in the yerrdown, yerrup arrays

        dyup = 2.303 * 10**phi * yerrup / 0.70**3 * h_cosmo**3
        zerrorup = 0.434 * dyup / 10**phi_plot
        dydown = 2.303 * 10**phi * yerrdown / 0.70**3 * h_cosmo**3
        zerrordown = 0.434 * dydown / 10**phi_plot

    ax1.errorbar(mass_plot, phi_plot, yerr = (zerrorup, zerrordown), label=allgalsobslabel, zorder=0, capsize=2)
    #  ax2.errorbar(mass, phi_new, yerr = (zerrorup, zerrordown), label=allgalsobslabel, zorder=0, capsize=2)




    colors = [['C5', 'C6', 'C7'],['C8', 'C9', 'C10']]

    # plot my simulations
    if not plot_obs_only:
        for i, simdir in enumerate(simdirs):
            allData, centralData, satelliteData, \
            orphanAllData, orphanSatData, filteredAllData, \
            filteredSatData = get_my_smf(simdir)

            # plot all without orphans
            mass, phi = allData
            ax1.plot(mass, phi, 
                    label=simname[simdir]+" centrals + satellites", 
                    c=colors[i][0])

            # plot all + orphan SMF
            mass, phi = orphanAllData
            ax1.plot(mass, phi, 
                    label=simname[simdir]+" centrals + satellites + orphans", 
                    c=colors[i][1])
            mass, phi = filteredAllData
            ax1.plot(mass, phi, ":", 
                    label=simname[simdir]+" centrals + satellites + filtered orphans", 
                    c=colors[i][2])

            # plot central galaxy SMF
            mass, phi = centralData
            ax2.plot(mass, phi, 
                    label=simname[simdir]+" centrals",
                    c=colors[i][0])

            # plot subhalo SMF
            mass, phi = satelliteData
            ax3.plot(mass, phi, 
                    label=simname[simdir]+" satellites", 
                    c=colors[i][0])

            # plot satellites and orphans
            mass, phi = orphanSatData
            ax3.plot(mass, phi, 
                    label=simname[simdir]+" satellites + orphans", 
                    c=colors[i][1])

            mass, phi = filteredSatData
            ax3.plot(mass, phi, ":", 
                    label=simname[simdir]+" satellites + filtered orphans",
                    c=colors[i][2])




    axtitle = r'$z \sim {0:4.2f}$'.format(redshift_chosen)

    for ax in fig.axes:
        ax.legend(title=axtitle, loc='lower left', prop=fontP)
        ax.set_ylim(-7.1,-0.2)

        if plot_obs_only:
            ax.set_ylabel(r"$\log_{10}\, \Phi(M_*) / (h^3 $Mpc$^{-3} $dex$^{-1})$")
            ax.set_xlabel(r"$\log_{10}\, M_* / (h^{-2} M_{\odot})$")
        else:
            ax.set_ylabel(r"$\log_{10}\, \Phi(M_*) / ($Mpc$^{-3} $dex$^{-1})$")
            ax.set_xlabel(r"$\log_{10}\, M_* / M_{\odot}$")
        #  ax.set_xlim(5,12.5)
        ax.grid()

    ax1.set_title("All Galaxies")
    ax2.set_title("Central Galaxies")
    ax3.set_title("Satellite Galaxies")

    plt.tight_layout()

    #  plt.show()
    if plot_obs_only:
        savefname = 'smf-obs-only.png'
        plt.savefig('smf-obs-only.png', dpi=300)
    else:
        savefname = 'smf-new-cosmo-with-satellites.png'
        plt.savefig(savefname, dpi=300)
    print("Saved", savefname)


    return









if __name__ == "__main__":
    main()
