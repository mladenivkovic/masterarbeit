#!/usr/bin/python3



import numpy as np
import warnings
import matplotlib.pyplot as plt
import Corrfunc.theory as ct
from Corrfunc.utils import convert_3d_counts_to_cf

import evaluation_module as em
from evaluation_module import g
from evaluation_module import ph



em.get_output_info()
em.read_info_files()
em.read_data()

rbin=np.logspace(-3, 1, 10)

np.random.seed(42)
N=g.galaxy_pos.shape[0]
boxsize=g.unit_l[g.ozi]/ph.Mpc
rand_N = 3*N
rand_X = np.random.uniform(0, boxsize, rand_N)
rand_Y = np.random.uniform(0, boxsize, rand_N)
rand_Z = np.random.uniform(0, boxsize, rand_N)
g.galaxy_pos *= boxsize
w = np.ones(N)
rand_w = np.ones(rand_N)


DD_counts=ct.DD(
    autocorr=1,
    nthreads=3,
    binfile=rbin,
    X1 = g.galaxy_pos[:,0],
    Y1 = g.galaxy_pos[:,1],
    Z1 = g.galaxy_pos[:,2],
    weights1 = w,
    #  weights1 = g.galaxy_masses,
    periodic=True,
    boxsize=boxsize
    )


DR_counts=ct.DD(
    autocorr=0,
    nthreads=3,
    binfile=rbin,
    X1 = g.galaxy_pos[:,0],
    Y1 = g.galaxy_pos[:,1],
    Z1 = g.galaxy_pos[:,2],
    #  weights1 = g.galaxy_masses,
    weights1 = w,
    periodic=True,
    boxsize=boxsize,
    X2 = rand_X,
    Y2 = rand_Y,
    Z2 = rand_Z,
    #  weights2 = g.galaxy_masses,
    weights2 = rand_w
    )

RR_counts=ct.DD(
    autocorr=1,
    nthreads=3,
    binfile=rbin,
    X1 = rand_X,
    Y1 = rand_Y,
    Z1 = rand_Z,
    #  weights1 = g.galaxy_masses,
    weights1 = rand_w,
    periodic=True,
    boxsize=boxsize
    )

for r in DD_counts: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10d} {4:10.6f}".
                         format(r['rmin'], r['rmax'], r['ravg'],
                         r['npairs'], r['weightavg']))

for r in DR_counts: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10d} {4:10.6f}".
                         format(r['rmin'], r['rmax'], r['ravg'],
                         r['npairs'], r['weightavg']))

for r in RR_counts: print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10d} {4:10.6f}".
                         format(r['rmin'], r['rmax'], r['ravg'],
                         r['npairs'], r['weightavg']))


cf = convert_3d_counts_to_cf(N, N, rand_N, rand_N,
        DD_counts, DR_counts,
        DR_counts, RR_counts)
for xi in cf: print("{0:10.6f}".format(xi))
