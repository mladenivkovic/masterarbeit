#!/usr/bin/python3

# Plots only particles belonging to defined halo. 

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

outputfilename = "plot_allhalo"



#==============
# Main routine
#==============
if __name__ == "__main__":

    # get particle data
    # halox: np.array
    halox, haloy, haloz = rf.get_all_halo_particles(noutput, nparticles) 

    fig, ax1 = rp.setup_fig()
    
    #plot halo first
    rp.plot_2D(halox, haloy, haloz, ax1, 'all', 0)


    rp.tweak_plot_2D(fig,ax1,'general')
    rp.save_fig(outputfilename,fig,workdir)
