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
halo=int(argv[2])
nparticles=int(argv[3])

outputfilename = "plot_halo_2D-"+str(halo)



#==============
# Main routine
#==============
if __name__ == "__main__":

    # get particle data
    # halox: np.array
    children, child_levels = rf.get_children_data(halo, noutput) 
    x_part, y_part, z_part, halox, haloy, haloz = rf.get_halo_particles(children, halo, noutput, nparticles)

    fig, ax1 = rp.setup_fig()
    
    #plot halo first
    rp.plot_2D(halox, haloy, haloz, ax1, 'halo', halo)

    #plot children
    for c,child in enumerate(children):
        if x_part[c].shape[0] > 10:
            rp.plot_2D(x_part[c], y_part[c],z_part[c], ax1, 'child', child)


    rp.tweak_plot_2D(fig,ax1,'general')
    rp.save_fig(outputfilename,fig,workdir)
