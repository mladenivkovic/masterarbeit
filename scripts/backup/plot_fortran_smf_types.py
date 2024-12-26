#!/usr/bin/env python3
#--------------------------------------------------
# Plots the stellar mass functions obtained from
# get_smf_types.f03.
# Needs smf-all.txt file as cmd line arg
# Assumes info_XXXXX.txt file is in the same dir
#--------------------------------------------------


import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()

import sys, os

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




#============================
def get_snapshot_data(fname):
#============================
    """
    Reads in info_XXXXX.txt file.
    """

    # find info_XXXXX.txt file
    dirname = os.path.dirname(fname)
    filelist = os.listdir(dirname)
    infofile = None
    for f in filelist:
        if 'info_' in f:
            infofile = f
            break

    if infofile is None:
        print("Didn't find infofile. Stopping.")
        quit(2)





    #-----------------------------
    # Read info files
    #-----------------------------

    try:
        ifile = open(os.path.join(dirname, infofile))
        for j in range(9):
            ifile.readline() # skip first 9 lines

        # get expansion factor
        aline = ifile.readline()
        astring, equal, aval = aline.partition("=")
        afloat = float(aval)
        a_exp = afloat

        for j in range(5):
            ifile.readline() # skip first 9 lines

        lline = ifile.readline()
        lstring, equal, lval=lline.partition('=')
        lfloat = float(lval)
        unit_l = lfloat/Mpc


    except IOError: # If file doesn't exist
        print("Didn't find any info data in ", dirname)

    redshift = 1.0/a_exp - 1

    return unit_l, redshift



def schechter(M, dlogm, phi_star, alpha, m_star):
    """
    """
    h = 0.6774
    m = M/10**(m_star)
    return phi_star*h**3*dlogm * m ** (alpha+1) * np.exp(-m)


def doubleschechter(M, dM):
    """
    For GAMA stuff
    https://arxiv.org/pdf/1111.5707.pdf
    """

    Mstar = 10**10.66
    phi1 = 3.96e-3
    phi2 = 0.79e-3
    alpha1 = -0.35
    alpha2 = -1.47
    MoM = M/Mstar

    res = np.exp(-MoM) * (phi1 * MoM** alpha1 + phi2 * MoM**alpha2) * dM/Mstar
    return res



#===================
def main():
#===================

    fname = sys.argv[-1]

    masses, central, satellite, orphan = np.loadtxt(fname,
            skiprows=1, dtype=np.float, unpack=True
            )

    unit_l, redshift = get_snapshot_data(fname)


    mass_center = 0.5*(masses[1:]+masses[:-1])
    dex_divide = masses[1:]-masses[:-1]
    masses = masses[:-1]
    central = central[:-1]
    satellite = satellite[:-1]
    orphan = orphan[:-1]

    # cut off at certain mass
    #  cutoff = mass_center > 7.8
    #  cutoff = np.logical_and(cutoff, mass_center<11.5)
    #  central = central[cutoff]
    #  satellite = satellite[cutoff]
    #  orphan = orphan[cutoff]
    #  mass_center = mass_center[cutoff]
    #  dex_divide = dex_divide[cutoff]
    #  masses = masses[cutoff]

    dex_divide_fit = np.zeros(mass_center.shape)
    dex_divide_fit = 0.1

    vol = (unit_l*(1+redshift))**3

    fig = plt.figure()
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)







    #  baldry = '/home/mivkov/UZH/Masterarbeit/masterarbeit/files/observational_data/GAMA_smf.txt'
    #  baldry_mass, baldry_phi = np.loadtxt(baldry, usecols=[0, 2], unpack=True)
    #  baldry_phi *= 1e-3





    #========================
    # Plot all
    #========================

    smf_all = central + satellite
    mask = smf_all > 0
    sma = smf_all[mask]/vol/dex_divide[mask]
    ax1.semilogy(mass_center[mask], sma, ls='--', label='excluding orphans')

    smf_all_orph = central + satellite + orphan
    mask = smf_all_orph > 0
    sma = smf_all_orph[mask]/vol/dex_divide[mask]
    ax1.semilogy(mass_center[mask], sma, ls='-', label='including orphans')


    masses_real = 10**masses
    dm = masses_real[1:]-masses_real[:-1]
    dlogm = masses[1:]-masses[:-1]


    #  fit = doubleschechter(masses_real[:-1], dm)/dex_divide[:-1]
    #  ax1.semilogy(mass_center[:-1], fit, label='fit GAMA doubleschechter')
    #  ax1.semilogy(baldry_mass, baldry_phi, label='baldry')


    fit = schechter(masses_real[:-1], dlogm, 0.01546, -1.164, 10.717) / dex_divide[:-1]
    ax1.semilogy(mass_center[:-1], fit, c='red', label='best schechter fit Yang+2009')


    #========================
    # Plot centrals
    #========================

    mask = central > 0
    smc = central[mask]/vol/dex_divide[mask]
    ax2.semilogy(mass_center[mask], smc, label='results')

    fit = schechter(masses_real[:-1], dlogm, 0.01084, -1.143, 10.758)/dex_divide[:-1]
    ax2.semilogy(mass_center[:-1], fit, c='red', label='best schechter fit Yang+2009')

    #  fit = doubleschechter(masses[:-1], 10**masses[1:]-10**masses[:-1])/dex_divide[:-1]
    #  ax2.semilogy(mass_center[:-1], fit, label='fit-test')


   #   mplot = np.loadtxt('/home/mivkov/Desktop/logM')
    #  phiplot = np.loadtxt('/home/mivkov/Desktop/centrals.txt')
    #  phierr = np.loadtxt('/home/mivkov/Desktop/centrals_errors.txt')
    #  h = 0.6774
    #  phiplot *= 1e-2 * h**3 / 0.1
    #  phierr *= 1e-2 * h**3
    #  #  ax2.errorbar(mplot, phiplot, yerr=phierr, label='data')
    #  ax2.semilogy(mplot, phiplot, label='data')


    #  ax2.semilogy(baldry_mass, baldry_phi, label='baldry')

    #========================
    # Plot satellites
    #========================

    mask = satellite > 0
    sms = satellite[mask]/vol/dex_divide[mask]
    ax3.semilogy(mass_center[mask], sms, label='excluding orphans')

    smso = satellite + orphan
    mask = smso > 0
    sms = smso[mask]/vol/dex_divide[mask]
    ax3.semilogy(mass_center[mask], sms, ls='--', label='including orphans')


    fit = schechter(masses_real[:-1], dlogm,  0.00692, -1.078, 10.483)/dex_divide[:-1]
    ax3.semilogy(mass_center[:-1], fit, c='red', label='best schechter fit Yang+2009')


    #  fit = doubleschechter(masses[:-1], 10**masses[1:]-10**masses[:-1])/dex_divide[:-1]
    #  ax3.semilogy(mass_center[:-1], fit, label='fit-test')



    ax1.set_title('All')
    ax2.set_title('Centrals')
    ax3.set_title('Satellites')

    ax1.set_xlim(8, 12)
    ax2.set_xlim(8, 12)
    ax3.set_xlim(8, 12)

    ax1.set_ylim(1e-6, 1e-1)
    ax2.set_ylim(1e-6, 1e-1)
    ax3.set_ylim(1e-6, 1e-1)

    ax1.legend()
    ax2.legend()
    ax3.legend()


    plt.show()



if __name__ == "__main__":
    main()
