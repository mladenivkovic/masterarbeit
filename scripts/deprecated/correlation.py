#!/usr/bin/python2


#=================================
# Calculate correlation function
#=================================



import numpy as np
import warnings
import matplotlib.pyplot as plt

import scipy.integrate

import evaluation_module as em
from evaluation_module import g
from evaluation_module import ph

nc = 100         # number of cells per dimension
mthresh = 0     # mass threshold [in M_Sol]


density_field = 0


#============================
def read_power_spectrum():
#============================
    """
    Reads in power spectrum from observations from file.
    """
    k_eff, k_low, k_high, P, dP = np.loadtxt('/home/mivkov/UZH/Masterarbeit/masterarbeit/observational_data/2dF_power_spectrum.txt',
        skiprows=4, unpack=True)

    k_high -= k_eff
    k_low = k_eff-k_low

    # change units from h Mpc-1 to Mpc-1

    k_high  /= g.h
    k_low   /= g.h
    k_eff   /= g.h
    P       *= g.h**3
    dP      *= g.h**3

    return k_eff, (k_low, k_high), P, dP






#====================
def main_v2():
#====================

    global g
    global density_field

    em.get_output_info()
    em.read_info_files()
    em.read_data()

    nc = 256
    l = g.unit_l[g.ozi]/ph.Mpc
    dk = 2*np.pi/l

    fig = plt.figure(figsize=(16,10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    axs = [ax1, ax2]

    cell_len = l/nc
    V = nc**3

    #  which = 'all'
    #  which = 'sub'
    which = 'main'

    density_field = em.get_density_field(nc, which=which)

    delta = density_field / np.mean(density_field) - 1.0

    delta_k = np.abs(np.fft.fftn(delta).round()) /V *(2*np.pi)**3
    Pk_field = delta_k**2




    #-----------------------------
    # Compute P(|k|)
    #-----------------------------

    # get 3d array of index integer distances to k = (0, 0, 0)
    dist = np.minimum(np.arange(nc), np.arange(nc,0,-1))
    dist *= dist
    dist_3d = np.sqrt(dist[:, None, None] + dist[:, None] + dist)

    distances_k, _ = np.unique(dist_3d, return_inverse=True)
    # add up how many counts you have, normalize by total number of counts to get mean
    # ravel: flatten ac array to 1D to correspond to _ array
    #  values_ac = np.bincount(_, weights=ac.ravel()) / np.bincount(_)
    Pk = np.bincount(_, weights=Pk_field.ravel()) / np.bincount(_)



    #-----------------------------
    # Compute xi(|r|)
    #-----------------------------

    xi_fft_field = np.fft.ifftn(Pk_field).real*(2*np.pi)**3

    dist_3d = np.sqrt(dist[:, None, None] + dist[:, None] + dist)
    distances_r, _ = np.unique(dist_3d, return_inverse=True)
    xi_fft = np.bincount(_, weights=xi_fft_field.ravel()) / np.bincount(_)

    dist_real = distances_r*cell_len
    k = distances_k*dk



    #  def to_integrate(r):
    #      return Pk * np.sin(k*r)*k/r
    #  xi_int = [scipy.integrate.trapz(to_integrate(r), k) for r in dist_real]

    def obs_correlation(r):
        #  https://arxiv.org/pdf/astro-ph/0301280.pdf
        return (r/(5.77/g.h))**(-1.80)

    ax1.semilogy(dist_real[xi_fft>0], xi_fft[xi_fft>0], label='xi_fft')
    #  ax1.semilogy(dist_real, xi_int, label='xi_int')
    ax1.semilogy(dist_real, obs_correlation(dist_real), label='obs')
    ax1.legend()

    k_obs, k_err, P_obs, P_err = read_power_spectrum()
    ax2.loglog(k, Pk)
    ax2.errorbar(k_obs, P_obs, yerr=P_err, xerr=k_err)

    plt.savefig('correlation-'+which+'.png')
    plt.show()


#====================
def main():
#====================

    global g
    global density_field

    em.get_output_info()
    em.read_info_files()
    em.read_data()

    nc =10
    l = g.unit_l[g.ozi]/ph.Mpc
    print "l=", l, "max", np.sqrt(3)*l/2
    cell_len = l/nc

    V = nc**3

    which = 'all'
    #  which = 'sub'
    #  which = 'main'

    density_field = em.get_density_field(nc, which=which)

    delta = density_field / np.mean(density_field) - 1.0
    print "Mean density:", np.mean(density_field)

    for i in range(nc):
        for j in range(nc):
            for k in range(nc):
                print "%12.4E" % delta[i,j,k],
            print
        print

    delta_k = np.abs(np.fft.rfftn(delta).round()) /V *(2*np.pi)**3
    Pk_field = delta_k**2




    #-----------------------------
    # Compute P(|k|)
    #-----------------------------

    # get 3d array of index integer distances to k = (0, 0, 0)
    dist = np.minimum(np.arange(nc), np.arange(nc,0,-1))
    dist_z = np.arange(nc//2+1)
    dist *= dist
    dist_z *= dist_z
    dist_3d = np.sqrt(dist[:, None, None] + dist[:, None] + dist_z)

    distances_k, _ = np.unique(dist_3d, return_inverse=True)
    # add up how many counts you have, normalize by total number of counts to get mean
    # ravel: flatten ac array to 1D to correspond to _ array
    # values_ac = np.bincount(_, weights=ac.ravel()) / np.bincount(_)
    Pk = np.bincount(_, weights=Pk_field.ravel()) / np.bincount(_)
    #  Pk_2, edges1 = np.histogram(dist_3d[:,:,1:-1], bins=distances_k, weights=Pk_field[:,:,1:-1])
    #  counts_2, edges2 = np.histogram(dist_3d[:,:,1:-1], bins=distances_k)


    # last bin is not included: it represents outermost limit.
    # I cut it off anyway.
    #  Pk_2[counts_2>0] /= counts_2[counts_2 >0]
    #  Pk[:-1] += Pk_2



    #-----------------------------
    # Compute xi(|r|)
    #-----------------------------

    xi_fft_field = np.fft.irfftn(Pk_field)*(2*np.pi)**3

    dist_3d = np.sqrt(dist[:, None, None] + dist[:, None] + dist)
    distances_r, _ = np.unique(dist_3d, return_inverse=True)
    xi_fft = np.bincount(_, weights=xi_fft_field.ravel()) / np.bincount(_)

    dist_real = distances_r*cell_len
    dk = 2*np.pi/l
    k = distances_k*dk





    #  def to_integrate(r):
    #      return Pk * np.sin(k*r)*k/r
    #  xi_int = [scipy.integrate.trapz(to_integrate(r), k) for r in dist_real]

    def obs_correlation(r):
        #  https://arxiv.org/pdf/astro-ph/0301280.pdf
        return (r/(5.77/g.h))**(-1.80)


    fig = plt.figure(figsize=(16,10))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.semilogy(dist_real[xi_fft>0], xi_fft[xi_fft>0], label='xi_fft')
    #  ax1.semilogy(dist_real, xi_int, label='xi_int')
    ax1.semilogy(dist_real, obs_correlation(dist_real), label='obs')
    ax1.legend()

    ax2.loglog(k, Pk)
    k_obs, k_err, P_obs, P_err = read_power_spectrum()
    ax2.errorbar(k_obs, P_obs, yerr=P_err, xerr=k_err)

    plt.savefig('correlation-'+which+'.png')
    plt.show()








#==================================
def main_simple():
#==================================


    nc = 128                # define how many cells your box has
    boxlen = 50.0           # define length of box
    Lambda = boxlen/4.0     # define an arbitrary wave length of a plane wave
    cellsize = boxlen/nc    # get size of a cell

    # create plane wave density field
    density_field = np.zeros((nc, nc, nc), dtype='float')
    for x in range(density_field.shape[0]):
        density_field[x,:,:] = np.cos(2*np.pi*x*cellsize/Lambda)

    # get overdensity field
    delta = density_field/np.mean(density_field) - 1

    # get P(k) field: explot fft of data that is only real, not complex
    delta_k = np.abs(np.fft.rfftn(delta).round()) #/(2.0*np.pi)**3
    Pk_field =  delta_k**2 #/ (2.0*np.pi)**3

    # get 3d array of index integer distances to k = (0, 0, 0)
    dist = np.minimum(np.arange(nc), np.arange(nc,0,-1))
    dist_z = np.arange(nc//2+1)
    dist *= dist
    dist_z *= dist_z
    dist_3d = np.sqrt(dist[:, None, None] + dist[:, None] + dist_z)

    # get unique distances and index which any distance stored in dist_3d
    # will have in "distances" array
    distances, _ = np.unique(dist_3d, return_inverse=True)

    # average P(kx, ky, kz) to P(|k|)
    Pk = np.bincount(_, weights=Pk_field.ravel())/np.bincount(_)

    # compute "phyical" values of k
    dk = 2*np.pi/boxlen
    k = distances*dk

    # plot results
    fig = plt.figure(figsize=(9,6))
    ax1 = fig.add_subplot(111)
    ax1.plot(k, Pk, label=r'$P(\mathbf{k})$')

    # plot expected peak:
    # k_peak = 2*pi/lambda, where we chose lambda for our planar wave earlier
    ax1.plot([2*np.pi/Lambda]*2, [Pk.min()-1, Pk.max()+1], label='expected peak')
    ax1.legend()
    plt.show()








#==================================
if __name__ == "__main__":
#==================================

    main()
