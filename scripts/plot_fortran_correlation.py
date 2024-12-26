#!/usr/bin/env python3

#====================================================================
# Plots the results of get_correlation.f03 scripts.
# This script is called by the get_correlation.sh script, but you can
# also call it manually:
#   $ plot_fortran_correlation.py
# it looks for the following files:
# correlation_all.txt, correlation_sub.txt
#====================================================================

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os.path import exists

redshift = 0
unit_l = 0
corr_all = 0
corr_sub = 0
outputlist = 0

Mpc = 3.086e24 # Mpc/cm
h = 0.6774



#============================
def get_snapshot_data():
#============================
    """
    Reads in correlation_*.txt and info_XXXXX.txt files.
    """

    from os import getcwd
    from os import listdir

    global redshift, unit_l, corr_all, corr_sub, outputlist

    workdir = getcwd()
    filelist = listdir(workdir)

    outputlist = []
    outdirtemp = 'output_'
    for filename in filelist:
        if outdirtemp in filename:
            # first check whether directory contains galaxy files
            txtfilelist = listdir(filename)
            if ('correlation_all.txt' in txtfilelist) and ('correlation_sub.txt' in txtfilelist):
                outputlist.append(filename)

    if len(outputlist)<1:
        print("I didn't find any output_XXXXX directories in current working directory.")
        print("Or none that contain an output_XXXXX/galaxies_XXXXX.txt00001 file.")
        print("Are you in the correct workdir?")
        quit()

    outputlist.sort()
    noutput = len(outputlist)
    print("First output containing correlation data: ", outputlist[0])

    a_exp = np.zeros(noutput, dtype='float')
    unit_l = np.zeros(noutput, dtype='float')

    #-----------------------------
    # Read info files
    #-----------------------------

    for i,out in enumerate(outputlist):
        fn = out+'/info_'+out[-5:]+'.txt'
        try:
            infofile = open(fn)
            for j in range(9):
                infofile.readline() # skip first 9 lines

            # get expansion factor
            aline = infofile.readline()
            astring, equal, aval = aline.partition("=")
            afloat = float(aval)
            a_exp[i] = afloat

            for j in range(5):
                infofile.readline() # skip first 9 lines

            lline = infofile.readline()
            lstring, equal, lval=lline.partition('=')
            lfloat = float(lval)
            unit_l[i] = lfloat/Mpc


        except IOError: # If file doesn't exist
            print("Didn't find any info data in ", srcdir)
            break

    redshift = 1.0/a_exp - 1




    #-------------------------------
    # Read smf.txt files
    #-------------------------------

    corr_all = [0 for i in outputlist]
    corr_sub = [0 for i in outputlist]
    for i,out in enumerate(outputlist):
        fn = out+'/correlation_all.txt'
        corr_all[i] = np.loadtxt(fn, dtype='float', skiprows=2,usecols=([0,1,2]))
        fn = out+'/correlation_sub.txt'
        corr_sub[i] = np.loadtxt(fn, dtype='float', skiprows=2,usecols=([0,1,2]))


    return



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
def plot_correlation():
#============================
    """
    plots correlations.
    """

    iz0 = np.argmin(np.absolute(redshift)) # find output closest to z=0
    r = corr_all[iz0][:,0]
    r = 0.5*(r[1:]+r[:-1])

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

    # Plot observational constraints
    ax.loglog(r, obs_correlation_1(r), '--', label='observational data 1')
    ax.loglog(r, obs_correlation_2(r), '--', label='observational data 2')


    corr_all_sum = np.zeros(r.shape[0],dtype='float')
    corr_all_counts_sum = np.zeros(r.shape[0],dtype='float')
    corr_sub_sum = np.zeros(r.shape[0],dtype='float')
    corr_sub_counts_sum = np.zeros(r.shape[0],dtype='float')
    for arr in corr_all:
        corr_all_sum += arr[1:,1]
        corr_all_counts_sum += arr[1:,2]
    for arr in corr_sub:
        corr_sub_sum += arr[1:,1]
        corr_sub_counts_sum += arr[1:,2]

    ztext = str(round(np.min(redshift),3))+r"$\leq$ z $\leq$ "+str(round(np.max(redshift),3))

    m = corr_all_counts_sum>0
    ax.loglog(r,corr_all_sum[m]/corr_all_counts_sum[m], label='including orphans,'+ztext)
    m = corr_sub_counts_sum>0
    ax.loglog(r,corr_sub_sum[m]/corr_sub_counts_sum[m], label='without orphans,'+ztext)


    ztext = r" z $\approx$ "+str(round(redshift[iz0],3))
    corr_all_z0 = corr_all[iz0][1:]
    corr_sub_z0 = corr_sub[iz0][1:]
    m = corr_all_z0[:,1]>0
    ax.loglog(r[m],corr_all_z0[:,1][m]/corr_all_z0[:,2][m], label='including orphans,'+ztext)
    m = corr_sub_z0[:,1]>0
    ax.loglog(r[m],corr_sub_z0[:,1][m]/corr_sub_z0[:,2][m], label='without orphans,'+ztext)



    ax.set_title("Correlation")
    ax.set_xlabel(r"$r$ $[Mpc (1+z)]$")
    ax.set_ylabel(r"$\xi(r)$ $[1]$")
    ax.legend()

    plt.savefig('correlation.png')




#============================
if __name__ == "__main__":
#============================
    get_snapshot_data()
    plot_correlation()


