#!/usr/bin/python3

#====================================================================
# Plots the results of get_correlation.f03 scripts.
# This script is called by the get_correlation.sh script, but you can
# also call it manually:
#   $ reportplot_correlations.py
# it looks for the following files:
# correlation_all.txt, correlation_sub.txt
# in the following directories:
#   galaxy-512/output_00041 and galaxy-512-100MPC/output_00067
#====================================================================

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif',
    'serif':['Computer Modern Roman'],
    'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)


srcdir = ''
rmax = 0



#=========================
def read_data(dname):
#=========================
    """
    Read cmdlineargs and data to plot
    """

    from sys import argv
    from os.path import exists

    #------------------------
    # Get case from cmdline
    #------------------------

    global srcdir,h

    xifile_all = dname+'/correlation_all.txt'
    xifile_sub = dname+'/correlation_sub.txt'

    if not exists(dname):
        print("Didn't find directory ", dname, "that you passed as argument. Aborting")
        quit()

    filelist=[xifile_all, xifile_sub]
    for f in filelist:
        if not exists(f):
            print("Didn't find file ", f)
            quit(2)


    r_sub,  xi_sub, xi_sub_counts = np.loadtxt(xifile_sub, dtype='float', unpack=True, skiprows=2)
    r_all,  xi_all, xi_all_counts = np.loadtxt(xifile_all, dtype='float', unpack=True, skiprows=2)

    data = [r_sub, xi_sub, xi_sub_counts, r_all, xi_all, xi_all_counts]


    dirnr = dname[-6:]
    infofname = dname+'/info'+dirnr+'.txt'
    infofile = open(infofname)
    for i in range(10):
        infofile.readline() # skip first 10 lines

    # get H0 factor
    hline = infofile.readline()
    hstring, equal, hval = hline.partition("=")
    hfloat = float(hval)
    h = hfloat/100.0
    infofile.close()


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
                if r[j] > rmax:
                    break
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
def plot_obs(x,y,label,ax):
#=============================================
    """
    Plots for Pk and xi
    """
    ax.loglog(x, y, label=label, zorder = 0)



#=============================================
def plot_wp(x,y,counts,label,ax,form='-'):
#=============================================
    """
    Plot projected correlation
    """
    mask = counts[1:]>0
    x_temp = 0.5*(x[1:]+x[:-1])
    ax.loglog(x_temp[mask],calc_wp(x_temp[mask], y[1:][mask]/counts[1:][mask]),form, label=label)






#==================================
def main():
#==================================

    #---------------------
    # Setup figure
    #---------------------

    fig1 = plt.figure(figsize=(6,6))
    fig2 = plt.figure(figsize=(6,6))
    ax1 = fig1.add_subplot(111)
    ax2 = fig2.add_subplot(111)

    global rmax



    #---------------------------------
    # Start plotting loop
    #---------------------------------

    dirs = [['galaxy-512-69MPC/output_00041', r'\texttt{G69} ', '--'],
            ['galaxy-512-100MPC/output_00041', r'\texttt{G100} ', ':']]

    for d in dirs:
        dname = d[0]
        addlabel = d[1]
        fmt = d[2]

        #-------------------------------
        # Read data
        #-------------------------------

        data = read_data(dname)

        # cut off ends because of too much scattering
        r_sub         = data[0][:-2]
        xi_sub        = data[1][:-2]
        xi_sub_counts = data[2][:-2]
        r_all         = data[3][:-2]
        xi_all        = data[4][:-2]
        xi_all_counts = data[5][:-2]

        rmaxtemp = max(r_sub[-1], r_all[-1])
        if rmax == 0:
            rmax = rmaxtemp
        else:
            rmax = min(rmax, rmaxtemp)


        #----------------------------------
        # Plot correlation
        #----------------------------------

        plot_stuff(r_sub, xi_sub, xi_sub_counts, addlabel+'excluding orphans', ax1, fmt)
        plot_stuff(r_all, xi_all, xi_all_counts, addlabel+'including orphans', ax1, fmt)


        #--------------------------------------
        # Plot projected correlation
        #--------------------------------------

        plot_wp(r_sub, xi_sub, xi_sub_counts, addlabel+'excluding orphans', ax2, fmt)
        plot_wp(r_all, xi_all, xi_all_counts, addlabel+'including orphans', ax2, fmt)





    #-----------------------------
    # Plot observations
    #-----------------------------

    plot_obs(r_all, obs_correlation_1(r_all), 'Zehavi et al 2004', ax1)
    plot_obs(r_all, obs_correlation_2(r_all), 'Li and White 2009', ax1)

    rp, wp, wperr= read_projected_corr()
    ax2.errorbar(rp, wp, yerr=wperr, label='Li and White 2009')



    #-----------------------------
    # Tweak plot, label axes
    #-----------------------------

    ax1.set_title("Correlation")
    ax1.set_xlabel(r"$r$ $[Mpc]$")
    ax1.set_ylabel(r"$\xi(r)$ $[1]$")
    ax1.legend()
    ax1.set_xlim(2.5e-2, 5e1)
    ax1.set_ylim(1e-4, 5e3)

    ax2.set_title("Projected Correlation")
    ax2.set_xlabel(r"$r_p$ $[Mpc]$ ")
    ax2.set_ylabel(r"$w_{p}(r_{p})$ $[Mpc]$")
    ax2.legend()


    figs = [(fig1, 'presentation-correlations'), (fig2, 'presentation-projected-correlations')]
    for f in figs:
        plt.figure(f[0].number)
        plt.tight_layout()

        plt.savefig(f[1]+".pdf", format='pdf')









#==================================
if __name__ == '__main__':
#==================================
    main()
