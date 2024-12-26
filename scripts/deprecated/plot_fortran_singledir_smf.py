#!/usr/bin/python3

#=========================================================
# Plots the results of singledir_eval_galaxies.f03
# scripts.
# This script is called by the singledir_eval_galaxies.sh
# script, but you can also call it manually:
# plot_fortran_correlation.py output_XXXXX
# it looks for "smf.txt" files
#=========================================================

import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from os.path import exists
from time import sleep









#------------------------
# Get case from cmdline
#------------------------


srcdir = argv[1]
smffile = srcdir+'/smf.txt'

if not exists(srcdir):
    print("Didn't find directory ", srcdir, "that you passed as argument. Aborting")
    quit()
if not exists(smffile):
    print("Didn't find file ", smffile)
    quit()


h = 0.704

#==============================
def get_observation_data():
#==============================
    """
    Read in data from saved file.
    """

    # https://arxiv.org/pdf/1111.5707.pdf 

    f = '/home/mivkov/UZH/Masterarbeit/masterarbeit/observational_data/GAMA_smf.txt'

    # log mass, bin width, number density, error, number in sample.
    
    logmass, binwidth, number_density, error, number_in_sample = np.loadtxt(f, dtype='float', skiprows=5, unpack=True)

    # number density is given in M_sol/1000 Mpc^-3; revert to Mpc
    number_density /= 1e3
    error /= 1e3

    return logmass, number_density, error
 




#-------------------------------
# Read and plot data
#-------------------------------

logmass, smf_all, smf_sub, smf_main = np.loadtxt(smffile, dtype='float', unpack=True, skiprows=1)

fig = plt.figure(figsize=(16,10))

ax1 = fig.add_subplot(111)
ax1.semilogy(logmass, smf_all, label='centrals, satellites & orphans')
ax1.semilogy(logmass, smf_sub, label='centrals & satellites')
ax1.semilogy(logmass, smf_main, label='centrals')
obs_logmass, obs_smf, obs_errors = get_observation_data()
ax1.errorbar(obs_logmass, obs_smf, yerr=obs_errors, label='observations')

ax1.set_title("Stellar Mass Function")
ax1.set_xlabel(r"$\log_{10}(M_{*})$")
ax1.set_ylabel(r"$\Psi$")
ax1.legend()

plt.savefig("smf.png")

