#!/usr/bin/python3

#------------------------------------------------------
# Plots the results of part2map.
# usage:
#   plot_part2map.py output_XXXXX
#   assumes output_XXXXX/part2map.dat file exists,
#   expects fortran binary, not ascii format.
#------------------------------------------------------

import numpy as np
import scipy.io
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from sys import argv


#================================
# Read cmdlineargs
#================================

if (len(argv)<2):
    print("I expect 1 argument: srcdir")
    quit(2)

srcdir = argv[1].strip()
if srcdir[-1] == '/':
    srcdir=srcdir[:-1]

mapfile = srcdir + '/part2map.dat'
if not os.path.exists(mapfile):
    print("I didn't find ", mapfile )
    quit(2)






#================================
# Read data
#================================

f = scipy.io.FortranFile(mapfile)

f.read_reals(dtype=np.float64) # t, dx, dy, dz
nx, ny = f.read_ints(dtype=np.int32)

data = f.read_reals(dtype=np.float32)
minval = data[data>0].min()
maxval = data[data>0].max()
data[data<minval] = minval*0.9 # cut off low end
data = data.reshape((nx,ny))

xmin, xmax = f.read_reals(dtype=np.float64)
ymin, ymax = f.read_reals(dtype=np.float64)



#=======================
# READ REDSHIFT
#=======================
fn = srcdir+'/info_'+srcdir[-5:]+'.txt'
try:
    infofile = open(fn)
    for j in range(9):
        infofile.readline() # skip first 9 lines

    # get expansion factor
    aline = infofile.readline()
    astring, equal, aval = aline.partition("=")
    afloat = float(aval)
    redshift = 1.0/afloat-1.0

    infofile.close()

    have_redshift = True

except IOError: # If file doesn't exist
    print("Didn't find any info data in ", srcdir, "won't be labelling redshift")
    have_redshift = False






#=============================
print("Creating Image")
#=============================


fig = plt.figure(figsize=(10,10), dpi=200)
ax=fig.add_subplot(111)

im = ax.imshow(data,
    interpolation='gaussian',
    cmap='afmhot',
    origin='lower',
    extent=(xmin,xmax,ymin,ymax),
    norm=LogNorm(),
    vmin=minval*0.9, vmax=maxval,
    )


# turn off axis
ax.set_axis_off()

# cut off margins
ax.xaxis.set_major_locator(plt.NullLocator())
ax.yaxis.set_major_locator(plt.NullLocator())

# annotate
ax.annotate(r'$z = {0:5.3f}$'.format(abs(redshift)), xy=(0.99,0.99), horizontalalignment='right', verticalalignment='top')





fig.tight_layout(pad=0.0,)

figname = 'particleplot_'+srcdir[-5:]+'.jpg'
print("saving figure ", figname)
plt.savefig(figname, form='jpg', quality=70, progressive=True, dpi=fig.dpi)
