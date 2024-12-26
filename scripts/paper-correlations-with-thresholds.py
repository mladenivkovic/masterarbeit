#!/usr/bin/env python3

#====================================================================
# Plots the results of get_correlation.f03 scripts.
# This script is called by the get_correlation.sh script, but you can
# also call it manually:
#   $ paper-correlations-with-thresholds.py
# it looks for the following files:
#   correlation_all-0.00E+00.txt, 
#   correlation_sub-0.00E+00.txt
#   correlation_all-1.00E+09.txt, 
#   correlation_sub-1.00E+09.txt
#   correlation_all-1.00E+10.txt, 
#   correlation_sub-1.00E+10.txt
#   correlation_all-1.00E+11.txt, 
#   correlation_sub-1.00E+11.txt
# in the following directories:
#   galaxy-512/output_00041 and galaxy-512-100MPC/output_00067
#====================================================================

import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
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

thresholds = ['0.00E+00', '1.00E+09', '1.00E+10', '1.00E+11']
threshold_labels = ["no threshold", r"$10^{9} M_\odot$", r"$10^{10} M_\odot$", r"$10^{11} M_\odot$", ]



#========================================
def read_data(dname, threshold=None):
#========================================
    """
    Read cmdlineargs and data to plot

    if threshold is none, read in all expected thresholds
    otherwise, read in just the one you provided as an index
    for the thresholds list
    """

    from sys import argv
    from os.path import exists

    #------------------------
    # Get case from cmdline
    #------------------------

    global srcdir,h

    if threshold is not None:
        thresholds_use = [thresholds[threshold]]
    else:
        thresholds_use = thresholds

    sub_data = []
    all_data = []

    for t in thresholds_use:
        xifile_all = dname+'/correlation_all-'+t+'.txt'
        xifile_sub = dname+'/correlation_sub-'+t+'.txt'

        if not exists(dname):
            print("Didn't find directory ", dname, "that you passed as argument. Aborting")
            quit()

        filelist=[xifile_all, xifile_sub]
        for f in filelist:
            if not exists(f):
                print("Didn't find file ", f)
                quit(2)


        print("reading", xifile_sub)
        r_sub,  xi_sub, xi_sub_counts = np.loadtxt(xifile_sub, dtype='float', unpack=True, skiprows=2)
        print("reading", xifile_all)
        r_all,  xi_all, xi_all_counts = np.loadtxt(xifile_all, dtype='float', unpack=True, skiprows=2)

        sub_data.append((r_sub, xi_sub, xi_sub_counts))
        all_data.append((r_all, xi_all, xi_all_counts))




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


    return sub_data, all_data





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






#=======================================
def read_projected_corr_observation():
#=======================================
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



#=========================================================
def plot_stuff(x,y,counts,label,ax,form={'linestyle':'-'}):
#=========================================================
    """
    Plots for Pk and xi
    """
    mask = counts[1:]>0
    x_temp = 0.5*(x[1:]+x[:-1])
    if isinstance(form, str):
        ax.loglog(x_temp[mask], y[1:][mask]/counts[1:][mask], form, label=label)
    else:
        ax.loglog(x_temp[mask], y[1:][mask]/counts[1:][mask], **form, label=label)
    return



#=============================================
def plot_obs(x,y,label,ax, fmt={'linestyle':'-'}):
#=============================================
    """
    Plots for Pk and xi
    """
    ax.loglog(x, y, label=label, zorder = 0, **fmt)



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
def plot_all_thresholds():
#==================================

    #---------------------
    # Setup figure
    #---------------------

    ncols = 2
    nrows = len(thresholds)

    fig = plt.figure(figsize=(4*ncols,4*nrows))
    for i in range(ncols*nrows):
        fig.add_subplot(nrows, ncols, i+1)



    global rmax



    #---------------------------------
    # Start plotting loop
    #---------------------------------

    #  dirs = [['galaxy-512-69MPC/output_00041', r'\texttt{G69} ', '--'],
    #          ['galaxy-512-100MPC/output_00041', r'\texttt{G100} ', ':']]
    dirs = [['galaxy-512-new-cosmo-69MPC/output_00067', r'\texttt{G69} ', '--'],
            ['galaxy-512-new-cosmo-100MPC/output_00067', r'\texttt{G100} ', ':']]

    for d in dirs:
        dname = d[0]
        addlabel = d[1]
        fmt = d[2]
        plot_observations = d == dirs[-1]

        #-------------------------------
        # Read data
        #-------------------------------

        sub_data, all_data = read_data(dname)


        #----------------------------------
        # Plot correlation
        #----------------------------------

        for i, t in enumerate(thresholds):

            sd = sub_data[i]
            ad = all_data[i]

            # cut off ends because of too much scattering
            r_sub         = sd[0][:-2]
            xi_sub        = sd[1][:-2]
            xi_sub_counts = sd[2][:-2]
            r_all         = ad[0][:-2]
            xi_all        = ad[1][:-2]
            xi_all_counts = ad[2][:-2]

            rmaxtemp = max(r_sub[-1], r_all[-1])
            if rmax == 0:
                rmax = rmaxtemp
            else:
                rmax = min(rmax, rmaxtemp)


            #--------------------------------
            # Plot my observations
            #--------------------------------

            ax = fig.axes[2*i] # even index; on the left side
            plot_stuff(r_sub, xi_sub, xi_sub_counts, addlabel+'excluding orphans', ax, fmt)
            plot_stuff(r_all, xi_all, xi_all_counts, addlabel+'including orphans', ax, fmt)

            #--------------------------------
            # Plot correlation observations
            #--------------------------------

            if plot_observations:
                plot_obs(r_all, obs_correlation_1(r_all), 'Zehavi et al 2004', ax)
                plot_obs(r_all, obs_correlation_2(r_all), 'Li and White 2009', ax)

            #--------------------------------------
            # Plot projected correlation
            #--------------------------------------

            ax = fig.axes[2*i+1] # odd index; on the right side
            plot_wp(r_sub, xi_sub, xi_sub_counts, addlabel+'excluding orphans', ax, fmt)
            plot_wp(r_all, xi_all, xi_all_counts, addlabel+'including orphans', ax, fmt)


            #----------------------------------------
            # Plot projected correlation observation
            #----------------------------------------
            if plot_observations:
                rp, wp, wperr= read_projected_corr_observation()
                ax.errorbar(rp, wp, yerr=wperr, label='Li and White 2009')



    #-----------------------------
    # Tweak plot, label axes
    #-----------------------------

    ax1 = fig.axes[0]
    ax1.set_title("Correlation")

    for i in range(nrows):
        ax = fig.axes[i*ncols]
        ax.set_xlabel(r"$r$ $[Mpc]$")
        ax.set_ylabel(r"$\xi(r)$ $[1]$")
        ax.legend(title=r"mass threshold: "+threshold_labels[i], fontsize='small')
        ax.set_xlim(9e-2, 5e1)
        ax.set_ylim(1e-4, 5e3)
        ax.grid()

    ax2 = fig.axes[1]
    ax2.set_title("Projected Correlation")
    for i in range(nrows):
        ax = fig.axes[i*ncols + 1]
        ax.set_xlabel(r"$r_p$ $[Mpc]$ ")
        ax.set_ylabel(r"$w_{p}(r_{p})$ $[Mpc]$")
        ax.legend(title=r"mass threshold: "+threshold_labels[i], fontsize='small')
        ax.grid()



    plt.tight_layout()

    plt.savefig("correlations-with-thresholds.pdf", format='pdf')






#==================================
def plot_one_threshold():
#==================================

    #---------------------
    # Setup figure
    #---------------------

    ncols = 1
    nrows = 1
    thresh = 1 # index in thresholds list to use!


    fig = plt.figure(figsize=(5*ncols,3.5*nrows))
    ax1 = fig.add_subplot(111)


    global rmax


    #---------------------------------
    # Start plotting loop
    #---------------------------------

    dirs = [['galaxy-512-new-cosmo-100MPC/output_00067', '', {'linestyle':'--', 'alpha':0.8}]]
    #  dirs = [['galaxy-512-new-cosmo-69MPC/output_00067', r'\texttt{G69} ', {'linestyle':'-', 'alpha':0.8}],
    #          ['galaxy-512-new-cosmo-100MPC/output_00067', r'\texttt{G100} ', {'linestyle':'--', 'alpha':0.8}]]

    for d in dirs:
        dname = d[0]
        addlabel = d[1]
        fmt = d[2]
        plot_observations = d == dirs[-1]

        #-------------------------------
        # Read data
        #-------------------------------

        sub_data, all_data = read_data(dname, threshold=thresh)


        #----------------------------------
        # Plot correlation
        #----------------------------------

        sd = sub_data[0]
        ad = all_data[0]

        # cut off ends because of too much scattering
        r_sub         = sd[0][:-2]
        xi_sub        = sd[1][:-2]
        xi_sub_counts = sd[2][:-2]
        r_all         = ad[0][:-2]
        xi_all        = ad[1][:-2]
        xi_all_counts = ad[2][:-2]

        rmaxtemp = max(r_sub[-1], r_all[-1])
        if rmax == 0:
            rmax = rmaxtemp
        else:
            rmax = min(rmax, rmaxtemp)


        #--------------------------------
        # Plot correlation observations
        #--------------------------------

        fmt2 = {"alpha":0.6}
        if plot_observations:
            r = np.logspace(-1.6, 1.52, 200)
            plot_obs(r, obs_correlation_1(r), 'Zehavi et al 2004', ax1, fmt2)
            plot_obs(r, obs_correlation_2(r), 'Li and White 2009', ax1, fmt2)

        #--------------------------------
        # Plot my observations
        #--------------------------------

        plot_stuff(r_sub, xi_sub, xi_sub_counts, addlabel+'excluding orphans', ax1, fmt)
        plot_stuff(r_all, xi_all, xi_all_counts, addlabel+'including orphans', ax1, fmt)




    #-----------------------------
    # Tweak plot, label axes
    #-----------------------------

    ax1 = fig.axes[0]
    ax1.set_title("Two-Point Correlation Function")
    ax1.legend(fontsize='small')

    ax1.set_xlabel(r"$r$ $[Mpc]$")
    ax1.set_ylabel(r"$\xi(r)$ $[1]$")
    ax1.set_xlim(3e-2, 7e1)
    ax1.set_ylim(1e-2, 5e4)
    ax1.grid()


    plt.tight_layout()

    plt.savefig("correlations-with-thresholds.pdf", format='pdf')









#==================================
if __name__ == '__main__':
#==================================
    plot_all_thresholds()
    #  plot_one_threshold()
