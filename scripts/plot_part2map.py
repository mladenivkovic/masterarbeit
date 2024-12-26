#!/usr/bin/python2
#------------------------------------------------------
# Plots the results of part2map.
# usage:
#   plot_part2map.py output_XXXXX
#   assumes output_XXXXX/part2map.dat file exists,
#   expects fortran binary, not ascii format.
#------------------------------------------------------

import numpy as np
import os.path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import fortranfile as ff
from sys import argv




#================================
# Read cmdlineargs
#================================

if (len(argv)<2):
    print "I expect 1 argument: srcdir"
    print "Optional argument: --plain to get only image, no axes, title or colorbar"
    quit(2)

srcdir = argv[1].strip()
if srcdir[-1] == '/':
    srcdir=srcdir[:-1]

mapfile = srcdir + '/part2map.dat'
if not os.path.exists(mapfile):
    print "I didn't find ", mapfile
    quit(2)

draw_labels = True

if len(argv)>=3:
    if argv[2] == '--plain' or argv[2] == '-p':
        draw_labels = False
        print "Drawing labels?", draw_labels
    else:
        print "Didn't recognize argument ", argv[2]
        print "Optional argument: --plain to get only image, no axes, title or colorbar"
        quit(2)





#=======================================
#  Reading FortranFile
#=======================================

f = ff.FortranFile(mapfile)

t, dx, dy, dz = f.readReals('d')
print "t=", t
nx, ny = f.readInts()


data = f.readReals()
#  print "SUM DATA:", sum(data), "mean:", np.mean(data), "min", data.min(), "max", data.max(), "nonzero-min", data[data>0].min()
#  print np.log10(data[data>0].min()), np.log10(data.max())

data[data<1e-8] = 1e-8 # cut off low end
data = data.reshape((nx,ny))

xmin, xmax = f.readReals('d')
ymin, ymax = f.readReals('d')






#=============================
print "Creating Image"
#=============================


fig = plt.figure(figsize=(10,10), dpi=300)
ax=fig.add_subplot(111)

im = ax.imshow(data,
    interpolation='gaussian',
    cmap='inferno',
    origin='lower',
    extent=(0,1,0,1),
    norm=LogNorm(),
    vmin=7e-9, vmax=3e-3
    )





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

except IOError: # If file doesn't exist
    print("Didn't find any info data in ", srcdir)


if draw_labels:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.15)
    fig.colorbar(im, cax=cax)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title("Particle Projection "+srcdir+' z={0:5.3f}'.format(redshift))
else:
    ax.set_title(r'$z=${0:5.3f}'.format(redshift))
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.set_yticklabels([])


fig.tight_layout()

figname = 'particleplot_'+srcdir[-5:]
if not draw_labels:
    figname += "_nolabels.pdf"
    print "saving figure ", figname
    plt.savefig(figname, form='pdf')
else:
    figname += ".png"
    print "saving figure ", figname
    plt.savefig(figname, form='png')


