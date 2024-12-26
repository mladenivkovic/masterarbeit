#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from friedmann import *

def a(z):
    return 1./ (1. + z)

def z(a, a0 = 1.):
    return a0/a - 1.

a_start = a(50)
a_end = 1.
nsnaps = 61

#  log_a_start = np.log(a_start)
#  log_a_end = np.log(a_end)
#
#  dloga = (log_a_end - log_a_start) / nsnaps
#  log_a = []
#
#  for i in range(nsnaps):
#      log_a.append(log_a_start + (i+1)*dloga)
#
#  log_a = np.array(log_a)
#  a_snaps = np.e ** log_a
#
#  print("{0:9} {1:9} {2:9}".format("Snapshot", "a", "z"))
#  for i in range(nsnaps):
#      print("{0:9d} {1:9.4f} {2:9.4f}".format(i+1, a_snaps[i], z(a_snaps[i])))

#  log_a_start = np.log10(a_start)
#  log_a_end = np.log10(a_end)
#  print(log_a_start, a_start, log_a_end, a_end)
#
#  log_a = np.logspace(log_a_start, log_a_end, nsnaps)
#  #  a_snaps = 10 ** log_a
#  a_snaps = log_a
#
#  print("{0:9} {1:9} {2:9}".format("Snapshot", "a", "z"))
#  for i in range(nsnaps):
#      print("{0:9d} {1:9.4f} {2:9.4f}".format(i+1, a_snaps[i], z(a_snaps[i])))
#
#  log_z_start = np.log(50)
#  log_z_end = np.log(1e-3)
#
#  dlogz = (log_z_end - log_z_start) / nsnaps
#  log_z = []
#
#  for i in range(nsnaps):
#      log_z.append(log_z_start + (i+1)*dlogz)
#
#  log_z = np.array(log_z)
#  z_snaps = np.e ** log_z
#
#  print("{0:9} {1:9} {2:9}".format("Snapshot", "a", "z"))
#  for i in range(nsnaps):
#      print("{0:9d} {1:9.4f} {2:9.4f}".format(i+1, a(z_snaps[i]), z_snaps[i]))

#  lin_a_start = a_start
#  lin_a_end = a_end
#
#  dlina = (lin_a_end - lin_a_start) / nsnaps
#  lin_a = []
#
#  for i in range(nsnaps):
#      lin_a.append(lin_a_start + (i+1)*dlina)
#
#  log_a = np.array(lin_a)
#  a_snaps = lin_a
#
#  print("{0:9} {1:9} {2:9}".format("Snapshot", "a", "z"))
#  for i in range(nsnaps):
#      print("{0:9d} {1:9.4f} {2:9.4f}".format(i+1, a_snaps[i], z(a_snaps[i])))


z_millenium = [127.0, 
79.997894, 
49.99959, 
30.000063, 
19.91569, 
18.243723, 
16.724525, 
15.343074, 
14.085914,
12.94078, 
11.89657, 
10.943864, 
10.073461, 
9.277915, 
8.549912, 
7.8832035, 
7.272188, 
6.7115865, 
6.1968336, 
5.723864, 
5.2888336, 
4.888449, 
4.5195556, 
4.1794686, 
3.8656828, 
3.575905, 
3.3080978, 
3.060419, 
2.8311827, 
2.6188614, 
2.422044, 
2.2394855, 
2.0700274, 
1.9126327, 
1.7663358, 
1.6302707, 
1.5036365, 
1.3857181, 
1.2758462, 
1.1734169, 
1.0778745, 
0.98870814, 
0.9054624, 
0.8276991, 
0.75503564, 
0.6871088, 
0.6235901, 
0.5641766, 
0.5085914, 
0.45657724, 
0.40789944, 
0.36234027, 
0.31970343, 
0.2798018, 
0.24246909, 
0.20754863, 
0.17489761, 
0.14438342, 
0.11588337, 
0.08928783, 
0.064493395, 
0.041403063, 
0.019932542, 
0.0] 

z_millenium = np.array(z_millenium)
z_millenium = z_millenium[2:]
nsnaps = z_millenium.shape[0]
a_millenium = a(z_millenium)

#  print("a millenium")
#  print(a_millenium)
#  times = get_times(a_millenium) + 13.8
#  print("times")
#  print(times)
#  times = times[::-1]
#  plt.figure()
#  plt.plot(times/times[-1], range(len(z_millenium)))
#  plt.show()


dloga = np.log10(a_millenium[-1]) - np.log10(a_millenium[-2])
additions = []
for i in range(3):
    additions.append(10**(np.log10(a_millenium[-1]) + dloga * (i+1) ))

full_a = list(a_millenium) + additions

for a in full_a:
    print("{0:.3f},".format(a), end="")
print()


