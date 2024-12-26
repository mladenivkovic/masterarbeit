#!/usr/bin/python3

#==================================================
# Plots particles belonging to a halo, different 
# clump levels are separated colour. 
#==================================================

from os import getcwd
from sys import argv #command line arguments

import ramsesplot_fetch as rf
import ramsesplot_plot2D as rp

#===========
# Set up
#===========
workdir= str(getcwd())
noutput=int(argv[1])
nparticles=int(argv[2])

outputfilename = "plot_levels"



#==============
# Main routine
#==============
if __name__ == "__main__":

    
    clump_level, halos = rf.get_level_data(noutput)
    x, y, z, hx, hy, hz = rf.get_all_level_particles(noutput, nparticles, clump_level, halos)


    fig, ax1 = rp.setup_fig()

    
    #plot halo first
    
    rp.plot_2D(hx, hy, hz, ax1, 'level-halo', 0)
    
    #plot subhaloes by level
    
    for level in range(len(x)-1,-1,-1):
        rp.plot_2D(x[level], y[level], z[level], ax1, 'level-subhalo', level)
    


    rp.tweak_plot_2D(fig,ax1,'general')
    rp.save_fig(outputfilename,fig,workdir)
