#!/usr/bin/python2


#=================================
# Calculate correlation function
#=================================



import numpy as np
import warnings
import matplotlib.pyplot as plt

import evaluation_module as em
from evaluation_module import g
from evaluation_module import ph

nc = 100        # number of cells per dimension
mthresh = 0     # mass threshold [in M_Sol]


density_field = 0



def read_power_spectrum():

    k_eff, k_low, k_high, P, dP = np.loadtxt('/home/mivkov/UZH/Masterarbeit/masterarbeit/observational_data/2dF_power_spectrum.txt',
        skiprows=4, unpack=True)

    k_high -= k_eff
    k_low = k_eff-k_low

    return k_eff, (k_low, k_high), P, dP






#====================
def main():
#====================

    global g
    global density_field

    em.get_output_info()
    em.read_info_files()
    em.read_data()

    which = 'main'
    density_field = em.get_density_field(nc, which=which)
    mean = np.mean(density_field)

    delta = (density_field / mean - 1.0 )

    ps = np.abs(np.fft.fftn(delta))
    ps *= ps

    ac = np.fft.ifftn(ps).real.round()
    ac /= ac[0, 0, 0]
    #  ps /=(2*np.pi)**3


    # Get distances for periodic box:
    # Increase integer(index) distance up until middle, then decrease
    dist = np.minimum(np.arange(nc), np.arange(nc, 0, -1))
    # gives [0, 1, 2, 3, 3, 1]
    # square it
    dist *= dist
    # create 3d map of distances to origin
    # for each dimension, add corresponding element in array for z to corresponding element for y for corresponding element for x
    # array elements are indices => integer distance from (0,0,0)
    dist_3d = np.sqrt(dist[:, None, None] + dist[:, None] + dist)
    # count which distances you have, and return an index array _ containing which index the distance in the distances array is for every element in dist_3d
    distances, _ = np.unique(dist_3d, return_inverse=True)
    # add up how many counts you have, normalize by total number of counts to get mean
    # ravel: flatten ac array to 1D to correspond to _ array
    values_ac = np.bincount(_, weights=ac.ravel()) / np.bincount(_)
    values_ps = np.bincount(_, weights=ps.ravel()) / np.bincount(_)
    

    l = g.unit_l[g.ozi]/ph.Mpc
    dist_real = distances[1:]*l/nc
    N = 1.0/nc**3
    k = 1.0/dist_real[::-1]

    def obs_correlation(r):
        return (r/(5.77/0.7))**(-1.80)

    fig = plt.figure(figsize=(16,10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    mask = values_ac>0
    mask[0] = False
    ax1.semilogy(distances[mask]*l/nc, values_ac[mask]/N/(2*np.pi)**3)
    ax1.semilogy(dist_real, obs_correlation(dist_real))
    
    ax2.loglog(k, values_ps[1:]*N)
    k_obs, k_err, P_obs, P_err = read_power_spectrum()
    ax2.errorbar(k_obs, P_obs, yerr=P_err, xerr=k_err)
    plt.show()
    
    plt.savefig('correlation-'+which+'.png')
    #  plt.show()
    #  plt.close()


    #  #  Test density field:
    #  fig = plt.figure()
    #  ax = fig.add_subplot(111)
    #
    #  ax.scatter(g.galaxy_pos[:,0], g.galaxy_pos[:,1], c='k', s=0.1, zorder=2)
    #  density_projection = np.sum(density_field, axis=2)
    #  # swap axes: imshow needs y axis to be first index
    #  density_projection = np.swapaxes(density_projection, 0, 1)
    #
    #  import matplotlib.colors as colors
    #  cellvolume = (g.unit_l[g.ozi]/ph.Mpc/nc)**3
    #  im = ax.imshow(density_projection,
    #      origin='lower',
    #      #  norm=colors.LogNorm(),
    #      cmap='autumn',
    #      extent=(0, 1, 0, 1),
    #      zorder=1)
    #
    #  plt.show()
    #


#==================================
if __name__ == "__main__":
#==================================

    main()
