#!/usr/bin/env python3

#============================================
# Compute and plot satellite fractions
# usage:
#   satellite_fraction.py output_XXXXX
#============================================


import numpy as np
from matplotlib import pyplot as plt


Mmin = 8
Mmax = 13
nbins = 50




def get_data():
    """
    Read in data
    return lists of masses for orphans, satellites and centrals
    """

    import sys
    import os

    dirname = sys.argv[1]
    if dirname[-1] ==  '/':
        dirname = dirname[:-1]

    dirnr = dirname[-5:]

    infofile = dirname+'/'+'info_'+dirname[-5:]+'.txt'
    f = open(infofile, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()

    ncpu = int(line[-1])
    f.close()


    halodata = [np.empty(1) for cpu in range(ncpu)]
    galaxydata = [np.empty((1, 5)) for cpu in range(ncpu)]

    for cpu in range(ncpu):
        hfile = dirname+'/halo_'+dirnr+'.txt'+str(cpu+1).zfill(5)
        halodata[cpu] = np.loadtxt(hfile, usecols=[0], skiprows=1)

        gfile = dirname+'/galaxies_'+dirnr+'.txt'+str(cpu+1).zfill(5)
        galaxydata[cpu] = np.loadtxt(gfile, usecols=[], skiprows=1)


    halos = np.concatenate(halodata)
    galaxies = np.concatenate(galaxydata, axis=0)

    which = np.zeros(galaxies.shape[0])

    # which:
    # 1 for halo, 2 for satellite, 0 for orphan
    nhalo = 0
    nsub = 0
    norph = 0

    for i in range(galaxies.shape[0]):
        if galaxies[i, 0] == 0:
            norph += 1
        else:
            if galaxies[i, 0] in halos:
                nhalo +=1
                which[i] = 1
            else:
                nsub += 1
                which[i] = 2

    print("nhalo:", nhalo, "nsub", nsub, "norph", norph)

    orphans = galaxies[:,1][which==0]
    halos = galaxies[:,1][which==1]
    subhalos = galaxies[:,1][which==2]


    return np.log10(orphans), np.log10(subhalos), np.log10(halos)



def main():

    orphans, satellites, centrals = get_data()


    hist_sub, bin_edges = np.histogram(satellites, range=(Mmin, Mmax), bins = nbins)
    hist_orph, bin_edges = np.histogram(orphans, range=(Mmin, Mmax), bins = nbins)
    hist_cent, bin_edges = np.histogram(centrals, range=(Mmin, Mmax), bins = nbins)

    bin_center = 0.5*(bin_edges[1:] + bin_edges[:-1])


    fract = hist_sub / (hist_sub + hist_cent)
    fract_orph = (hist_sub+hist_orph)/(hist_sub+hist_orph+hist_cent)


    fig = plt.figure()
    ax = fig.add_subplot(111)

    #  ax.plot(bin_center, hist_sub, label='sub')
    #  ax.plot(bin_center, hist_orph, label='orphan')
    #  ax.plot(bin_center, hist_cent, label='halo')

    ax.plot(bin_center, fract, label='without orphans')
    ax.plot(bin_center, fract_orph, label='with orphans')



    ax.legend()

    plt.show()





if __name__ == '__main__':
    main()
