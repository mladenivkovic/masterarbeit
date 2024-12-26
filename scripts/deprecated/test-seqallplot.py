#!/usr/bin/python3

#==================================================
# Create plot that is part of a sequence for 
# multiple output_* directories
#==================================================

from os import getcwd
from sys import argv #command line arguments

import ramsesplot_fetch as rf
import ramsesplot_plot2D as rp
import testplot_fetch as tf
import testplot_plot as tp

import numpy as np


#===========
# Set up
#===========
workdir= str(getcwd())
noutput=int(argv[1])
nparticles=int(argv[2])

temp_name=argv[4] # filename of first output_*/clump_* 

outputfilename = "png_seq/all_seqplot_"+temp_name[7:12]




#==============
# Main routine
#==============
if __name__ == "__main__":
    
    # get halo
    halos, subhalos = tf.get_all_clump_data(4, noutput)

    halox, haloy, haloz, subhalox, subhaloy, subhaloz = tf.get_all_clump_particles(noutput, nparticles, halos, subhalos)

    fig, ax1, ax2, ax3 = rp.setup_3fig()
    
    #plot halo first
    colorindex = 0
    for i, halo in enumerate(halos):
        tp.plot_32D(np.asarray(halox[i]), np.asarray(haloy[i]), np.asarray(haloz[i]), ax1, ax2, ax3, 'seq-halo', halo, colorindex)
        colorindex += 1

    #plot subhalos
    for i, subhalo in enumerate(subhalos):
        if len(subhalox[i]) > 0:
            colorindex += 1
            tp.plot_32D(np.asarray(subhalox[i]), np.asarray(subhaloy[i]), np.asarray(subhaloz[i]), ax1, ax2, ax3, 'seq-child', subhalo, colorindex)


    rp.tweak_3plot_2D(fig, ax1, ax2, ax3, 'general')
    rp.save_fig(outputfilename, fig, workdir)

