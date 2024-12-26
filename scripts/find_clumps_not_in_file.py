#!/usr/bin/env python3

#---------------------------------------------------------
# find_clumps_not_in_file.py file1 file2
# find clumps in file2 that are not in file1
# 
# This part is to identify stong fluctuations
# between different simulations.
# Extract descendants and progenitors, which you then
# should copypaste into run_extreme_growth_plots.sh
# file to run the comparisons. The input files are 
# intended to be stdout output from evaltrees.py
#---------------------------------------------------------

import numpy as np
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)

snap1 = data1[:,1]
snap2 = data2[:,1]

snapmin = int(snap2.min())
snapmax = int(snap2.max())

tol = 1e-2
nmatch = 0

unmatched_desc = []
unmatched_prog = []
unmatched_desc_snap = []
unmatched_prog_snap = []
unmatched_desc_mass = []
unmatched_prog_mass = []

for s in range(snapmin, snapmax+1):
    mask1 = snap1 == s
    mask2 = snap2 == s

    desc1 = data1[mask1,0]
    dmass1 = data1[mask1,2]
    xd1 = data1[mask1,3]
    yd1 = data1[mask1,4]
    zd1 = data1[mask1,5]
    prog1 = data1[mask1,6]
    psnap1 = data1[mask1,7]
    pmass1 = data1[mask1, 8]
    xp1 = data1[mask1,9]
    yp1 = data1[mask1,10]
    zp1 = data1[mask1,11]
    dnparts1 = data1[mask1,12]

    desc2 = data2[mask2,0]
    dmass2 = data2[mask2,2]
    xd2 = data2[mask2,3]
    yd2 = data2[mask2,4]
    zd2 = data2[mask2,5]
    prog2 = data2[mask2,6]
    psnap2 = data2[mask2,7]
    pmass2 = data2[mask2, 8]
    xp2 = data2[mask2,9]
    yp2 = data2[mask2,10]
    zp2 = data2[mask2,11]
    dnparts2 = data2[mask2,12]

    has_match1 = np.zeros(desc1.shape[0])
    has_match2 = np.zeros(desc2.shape[0])

    for i in range(desc2.shape[0]):
        
        for j in range(desc1.shape[0]):
            if not has_match1[j]:
                dxd = (xd1[j] - xd2[i])**2 + (yd1[j] - yd2[i])**2 + (zd1[j] - zd2[i])**2
                dxd = np.sqrt(dxd)
                dxp = (xp1[j] - xp2[i])**2 + (yp1[j] - yp2[i])**2 + (zp1[j] - zp2[i])**2
                dxp = np.sqrt(dxp)

                is_match = abs(dxd) <= tol
                is_match = is_match and abs(dxp) <= tol
                is_match = is_match and abs(1. - dmass2[i]/dmass1[j]) <= tol

                #  print(abs(dxd), abs(dxp), abs(1. - dmass2[i]/dmass1[j]))
                if is_match:
                    has_match1[j] = 1
                    has_match2[i] = 1
                    #  print("Found match file1", desc1[j], "file2", desc2[i])
                    nmatch += 1
                    break

    unmatched = np.logical_not(has_match2)
    #  print(desc2)
    #  print(desc2[unmatched])
    #  print(unmatched)
    #  quit()
    unmatched_desc.append(desc2[unmatched])
    unmatched_prog.append(prog2[unmatched])
    unmatched_desc_snap.append(snap2[mask2][unmatched])
    unmatched_prog_snap.append(psnap2[unmatched])
    unmatched_desc_mass.append(dmass2[unmatched])
    unmatched_prog_mass.append(pmass2[unmatched])

print("found", nmatch, "/", data2.shape[0], "matches")

unmatched_desc = np.concatenate(unmatched_desc).astype(int)
unmatched_prog = np.concatenate(unmatched_prog).astype(int)
unmatched_desc_snap = np.concatenate(unmatched_desc_snap).astype(int)
unmatched_prog_snap = np.concatenate(unmatched_prog_snap).astype(int)
unmatched_desc_mass = np.concatenate(unmatched_desc_mass)
unmatched_prog_mass = np.concatenate(unmatched_prog_mass)

progsort = np.argsort(unmatched_prog_mass)[::-1]
descsort = np.argsort(unmatched_desc_mass)[::-1]

for i in range(20):
    dind = descsort[i]
    pind = progsort[i]
    print(unmatched_desc[dind], unmatched_desc[pind], end=" ")
print()

for i in range(20):
    dind = descsort[i]
    pind = progsort[i]
    print(unmatched_desc_snap[dind], unmatched_desc_snap[pind], end=" ")
print()

for i in range(20):
    dind = descsort[i]
    pind = progsort[i]
    print(unmatched_prog[dind], unmatched_prog[pind], end=" ")
print()

for i in range(20):
    dind = descsort[i]
    pind = progsort[i]
    print(unmatched_prog_snap[dind], unmatched_prog_snap[pind], end=" ")
print()

