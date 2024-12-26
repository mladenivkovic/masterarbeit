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


#===========
# Set up
#===========
workdir= str(getcwd())
noutput=int(argv[1])
nparticles=int(argv[2])

temp_name=argv[4] # filename of first output_*/clump_* 

outputfilename = "png_seq/seqplot_"+temp_name[7:12]




#==============
# Main routine
#==============
if __name__ == "__main__":
    
    # get halo
    halo=tf.find_halo(4, noutput)

    children, child_levels = rf.get_children_data(halo, noutput) 
    x_part, y_part, z_part, halox, haloy, haloz = rf.get_halo_particles(children, halo, noutput, nparticles)


    fig, ax1, ax2, ax3 = rp.setup_3fig()
    
    #plot halo first
    colorindex = 0
    tp.plot_32D(halox, haloy, haloz, ax1, ax2, ax3, 'seq-halo', halo, colorindex)

    #plot children
    for c,child in enumerate(children):
        if x_part[c].shape[0] > 10:
            colorindex += 1
            tp.plot_32D(x_part[c], y_part[c], z_part[c], ax1, ax2, ax3, 'seq-child', child, colorindex)


    rp.tweak_3plot_2D(fig, ax1, ax2, ax3, 'general')
    rp.save_fig(outputfilename, fig, workdir)

