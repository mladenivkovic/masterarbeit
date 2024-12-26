#!/usr/bin/python2


#====================================================
#
# Plot halo mass function and stellar mass function
# at z = 0
#
#=====================================================


import numpy as np
import warnings
import matplotlib.pyplot as plt

import evaluation_module as em
from evaluation_module import g
from evaluation_module import ph








#============================
def get_smf():
#============================
    """
    Calculate the stellar mass function for given output
    oind:  index of output
    """

    bins = [ 6.0, 6.5,
        7.0, 7.2, 7.4, 7.6, 7.8,
        8.0, 8.2, 8.4, 8.6, 8.8,
        9.0, 9.2, 9.4, 9.6, 9.8,
        10.0, 10.2, 10.4, 10.6, 10.8,
        11.0, 11.2, 11.4, 11.6, 11.8,
        12.0
    ]


    galaxy_masses = np.log10(g.galaxy_masses)

    # Get stellar masses of main halos
    sm = np.zeros(galaxy_masses.shape[0], dtype='float')
    sm_sub = np.zeros(galaxy_masses.shape[0], dtype='float')
    sm_sub_orph = np.zeros(galaxy_masses.shape[0], dtype='float')

    i1 = 0
    i2 = 0
    i3 = 0
    for j,c in enumerate(g.galaxy_clumps):
        if c > 0:  # exclude orphans
            if c in g.halos: # exclude subhaloes
                sm[i1] = galaxy_masses[j]
                i1+= 1
            sm_sub[i2] = galaxy_masses[j]
            i2 += 1
        sm_sub_orph[i3] = galaxy_masses[j]
        i3 += 1


    sm = sm[:i1+1]
    sm_sub = sm_sub[:i2+1]
    sm_sub_orph = sm_sub_orph[:i3+1]


    smf, binedges = np.histogram(sm, bins=bins)
    smf_sub, binedges = np.histogram(sm_sub, bins=bins)
    smf_sub_orph, binedges = np.histogram(sm_sub_orph, bins=bins)

    #  smf = smf / (10.0**binedges[1:]-10.0**binedges[:-1])
    smf = smf / (binedges[1:]-binedges[:-1])
    smf_sub = smf_sub / (binedges[1:]-binedges[:-1])
    smf_sub_orph = smf_sub_orph / (binedges[1:]-binedges[:-1])


    vol = (g.unit_l[g.ozi]/10/ph.Mpc)**3 # observational number density is per 10^3 Mpc
    smf = smf/vol
    smf_sub = smf_sub/vol
    smf_sub_orph = smf_sub_orph/vol
    masses = (binedges[1:]+binedges[:-1])/2


    return smf, smf_sub, smf_sub_orph, masses







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

    return number_density, error







#====================
def main():
#====================

    global g

    em.get_output_info()
    em.read_info_files()
    em.read_data(read_pos=False)


    smf, smf_sub, smf_sub_orph, masses = get_smf()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    obs_smf, obs_errors = get_observation_data()
    ax.errorbar(masses, obs_smf, yerr=obs_errors, label='observations')
    ax.semilogy(masses, smf, label='centrals only')
    ax.semilogy(masses, smf_sub, label='subhalos')
    ax.semilogy(masses, smf_sub_orph, label='subhalos + orphans')
    ax.legend()




    plt.show()



#==================================
if __name__ == "__main__":
#==================================

    main()
