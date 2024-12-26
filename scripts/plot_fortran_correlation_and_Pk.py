#!/usr/bin/python3

#====================================================================
# Plots the results of get_correlation.f03 scripts.
# This script is called by the get_correlation.sh script, but you can
# also call it manually:
#   $ plot_fortran_correlation_and_Pk.py output_XXXXX
# it looks for the following files:
# Pk_all.txt, Pk_sub.txt, correlation_all.txt, correlation_sub.txt
#====================================================================

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                                'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)


h = 0.704
#  h = 0.6774

srcdir = ''



#=========================
def read_data():
#=========================
    """
    Read cmdlineargs and data to plot
    """

    from sys import argv
    from os.path import exists

    #------------------------
    # Get case from cmdline
    #------------------------

    global srcdir

    try:
        srcdir = argv[1]
    except IndexError:
        print("I need the output_XXXXX directory as argument.")
        quit()



    pkfile_all = srcdir+'/Pk_all.txt'
    pkfile_sub = srcdir+'/Pk_sub.txt'
    xifile_all = srcdir+'/correlation_all.txt'
    xifile_sub = srcdir+'/correlation_sub.txt'

    if not exists(srcdir):
        print("Didn't find directory ", srcdir, "that you passed as argument. Aborting")
        quit()

    filelist=[pkfile_all, pkfile_sub, xifile_all, xifile_sub]
    for f in filelist:
        if not exists(f):
            print("Didn't find file ", f)
            quit(2)


    k_sub,  Pk_sub, Pk_sub_counts = np.loadtxt(pkfile_sub, dtype='float', unpack=True, skiprows=2)
    k_all,  Pk_all, Pk_all_counts = np.loadtxt(pkfile_all, dtype='float', unpack=True, skiprows=2)
    r_sub,  xi_sub, xi_sub_counts = np.loadtxt(xifile_sub, dtype='float', unpack=True, skiprows=2)
    r_all,  xi_all, xi_all_counts = np.loadtxt(xifile_all, dtype='float', unpack=True, skiprows=2)
    #  rp_sub, wp_sub, wp_sub_counts = np.loadtxt(wpfile_sub, dtype='float', unpack=True, skiprows=2)
    #  rp_all, wp_all, wp_all_counts = np.loadtxt(wpfile_all, dtype='float', unpack=True, skiprows=2)


    data = [k_sub, Pk_sub, Pk_sub_counts, k_all, Pk_all, Pk_all_counts,
        r_sub, xi_sub, xi_sub_counts, r_all, xi_all, xi_all_counts]

    return data





#=========================
def obs_correlation_1(r):
#=========================
    #  https://arxiv.org/pdf/astro-ph/0301280.pdf
    return (r/(5.77/h))**(-1.80)
#=========================
def obs_correlation_2(r):
#=========================
    #  https://arxiv.org/pdf/0901.0706.pdf
    return (r/6.1*h)**(-1.84)




#============================
def read_power_spectrum():
#============================
    """
    Reads in power spectrum from observations from file.
    """
    obsfile='/home/mivkov/UZH/Masterarbeit/masterarbeit/files/observational_data/2dF_power_spectrum.txt'
    k_eff, k_low, k_high, P, dP = np.loadtxt(obsfile, skiprows=4, unpack=True)

    dk = k_high - k_low
    k_high -= k_eff
    k_low = k_eff-k_low


    # change units from h Mpc-1 to Mpc-1

    k_high  *= h
    k_low   *= h
    k_eff   *= h
    dk      *= h
    P       /= h**3
    dP      /= h**3

    #  return k_eff[1:], (k_low[1:], k_high[1:]), P[1:]/np.log10(dk[1:]/dk[:-1]), dP[1:]/np.log10(dk[1:]/dk[:-1])
    return k_eff, (k_low, k_high), P, dP



#============================
def read_projected_corr():
#============================
    """
    Reads in mass weighted projected correlation from Li&White 2009
    """
    obsfile='/home/mivkov/UZH/Masterarbeit/masterarbeit/files/observational_data/projected_corr_LiWhite.txt'
    rp, wp, wperr = np.loadtxt(obsfile,skiprows=1,unpack=True)


    # change units from h-1 Mpc to Mpc

    rp /= h
    wp /= h
    wperr /= h

    return rp, wp, wperr



#=============================
def calc_wp(r, xi):
#=============================
    """
    Calculate projected corr by integrating correlation function
    """

    def get_integrand(j):
        return r[j]*xi[j]/np.sqrt(r[j]**2-rp**2)

    wp = np.zeros(r.shape[0],dtype='float')

    for i, rp in enumerate(r[:-1]):
        # special case for first
        dr = r[i+1]-r[i]
        wp[i] += 2*get_integrand(i+1)*dr

        if i+2 < r.shape[0]:
            for j in range(i+2, r.shape[0]):
                dr = r[j]-r[j-1]
                # factor 2 and 0.5 for mean [ 0.5 * (int[j]+int[j-1]) * dr ] cancel out
                wp[i] += (get_integrand(j)+get_integrand(j-1))*dr

    return wp



#=============================================
def plot_stuff(x,y,counts,label,ax,form='-'):
#=============================================
    """
    Plots for Pk and xi
    """
    mask = counts[1:]>0
    x_temp = 0.5*(x[1:]+x[:-1])
    ax.loglog(x_temp[mask], y[1:][mask]/counts[1:][mask], form, label=label)


#=============================================
def plot_wp(x,y,counts,label,ax,form='-'):
#=============================================
    """
    Plot projected correlation
    """
    mask = counts>0
    ax.loglog(x[mask],calc_wp(x[mask], y[mask]/counts[mask]),form, label=label)






#==================================
def main():
#==================================


    #-------------------------------
    # Read data
    #-------------------------------

    data = read_data()

    k_sub         = data[0]
    Pk_sub        = data[1]
    Pk_sub_counts = data[2]
    k_all         = data[3]
    Pk_all        = data[4]
    Pk_all_counts = data[5]
    r_sub         = data[6]
    xi_sub        = data[7]
    xi_sub_counts = data[8]
    r_all         = data[9]
    xi_all        = data[10]
    xi_all_counts = data[11]



    fig = plt.figure(figsize=(16,6))

    #-------------------------------
    # Plot power spectrum
    #-------------------------------

    ax1 = fig.add_subplot(131)

    plot_stuff(k_sub, Pk_sub, Pk_sub_counts, 'excluding orphans', ax1)
    plot_stuff(k_all, Pk_all, Pk_all_counts, 'including orphans', ax1)

    # Observational Data
    k_obs, k_err, P_obs, P_err = read_power_spectrum()
    ax1.errorbar(k_obs, P_obs, yerr=P_err, xerr=k_err, label='observational data')


    #----------------------------------
    # Plot correlation
    #----------------------------------

    ax2 = fig.add_subplot(132)

    plot_stuff(r_sub, xi_sub, xi_sub_counts, 'excluding orphans', ax2)
    plot_stuff(r_all, xi_all, xi_all_counts, 'including orphans', ax2)
    plot_stuff(r_all[1:], obs_correlation_1(r_all[1:]), np.ones(r_all.shape[0]-1), 'Zehavi et al 2004', ax2)
    plot_stuff(r_all[1:], obs_correlation_2(r_all[1:]), np.ones(r_all.shape[0]-1), 'Li and White 2009', ax2)




    #--------------------------------------
    # Plot projected correlation
    #--------------------------------------

    ax3 = fig.add_subplot(133)
    x = np.linspace(1,100,1000)

    plot_wp(r_sub, xi_sub, xi_sub_counts, 'excluding orphans', ax3)
    plot_wp(r_all, xi_all, xi_all_counts, 'including orphans', ax3)

    rp, wp, wperr= read_projected_corr()
    ax3.errorbar(rp, wp, yerr=wperr, label='Li and White 2009')



    #-----------------------------
    # Tweak plot, label axes
    #-----------------------------

    ax1.set_title("Power Spectrum")
    ax1.set_xlabel(r"$k$ $[Mpc^{-1}]$")
    ax1.set_ylabel(r"$P(k)$ $[Mpc^{3}]$")
    ax1.legend()

    ax2.set_title("Correlation")
    ax2.set_xlabel(r"$r$ $[Mpc]$")
    ax2.set_ylabel(r"$\xi(r)$ $[1]$")
    ax2.legend()

    ax3.set_title("Projected Correlation")
    ax3.set_xlabel(r"$r_p$ $[Mpc]$ ")
    ax3.set_ylabel(r"$w_{p}(r_{p})$ $[Mpc]$")
    ax3.legend()



    plt.tight_layout()

    plt.savefig("correlation_"+srcdir[-5:]+".pdf", format='pdf')









#==================================
if __name__ == '__main__':
#==================================
    main()
