#!/usr/bin/python3

#==================================================
# Plot halos by levels, mark most bound particles
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

outputfilename = "mostboundplot_levels"+temp_name[7:12]



#==============
# Main routine
#==============
if __name__ == "__main__":

    # get levels 
    clump_level, halos = rf.get_level_data(noutput)
    px, py, pz, hx, hy, hz, mbx, mby, mbz, mbhx, mbhy, mbhz = tf.get_all_level_particles(noutput, nparticles, clump_level, halos)


    fig, ax1, ax2, ax3 = rp.setup_3fig()

    colorindex = 0
    #plot halo first
    tp.plot_32D(hx, hy, hz, ax1, ax2, ax3, 'level-halo', 0, colorindex)
    tp.plot_32D(mbhx, mbhy, mbhz, ax1, ax2, ax3, 'level-halo-mb', 0, colorindex)

    #plot other levels
    for level in range(len(px)-1,-1,-1):
        colorindex += 1
        tp.plot_32D(px[level], py[level], pz[level], ax1, ax2, ax3, 'level-subhalo', level, colorindex)
        tp.plot_32D(mbx[level], mby[level], mbz[level], ax1, ax2, ax3, 'level-subhalo-mb', level, colorindex)

   

    rp.tweak_3plot_2D(fig, ax1, ax2, ax3, 'general')
    rp.save_fig(outputfilename, fig, workdir)
