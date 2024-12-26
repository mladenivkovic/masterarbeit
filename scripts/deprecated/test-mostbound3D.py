#!/usr/bin/python3

#==================================================
# Plots only particles belonging to defined halo. 
# marks most bound particles 
#==================================================

from os import getcwd
from sys import argv #command line arguments

import ramsesplot_fetch as rf
import ramsesplot_plot3D as rp
import testplot_fetch as tf
import testplot_plot as tp

#===========
# Set up
#===========
workdir= str(getcwd())
noutput=int(argv[1])
halo=int(argv[2])
nparticles=int(argv[3])

outputfilename = "mostboundplot_halo_"+str(halo)



#==============
# Main routine
#==============
if __name__ == "__main__":

    # Read in halo data: Find children

    # get clump data
    children, child_levels = rf.get_children_data(halo,noutput)

    # get particle data
    # x_part: list of np.arrays
    # halox: np.array
    x_part, y_part, z_part, most_unb_x, most_unb_y, most_unb_z, halox, haloy, haloz, halo_mbx, halo_mby, halo_mbz = tf.get_halo_particles(children, halo, noutput, nparticles)
    

    fig, ax1 = rp.setup_fig()
    
    colorindex = 0
    #plot halo first
    tp.plot_3D(halox, haloy, haloz, ax1, 'halo-transp', halo, colorindex)
    tp.plot_3D(halo_mbx, halo_mby, halo_mbz, ax1, 'halo-mb', halo, colorindex)

    #plot children
    for c,child in enumerate(children):
        if len(x_part[c]) > 10:
            colorindex += 1
            tp.plot_3D(x_part[c], y_part[c],z_part[c], ax1, 'child-transp', child, colorindex)
            tp.plot_3D(most_unb_x[c], most_unb_y[c], most_unb_z[c], ax1, 'child-mb', child, colorindex)


    rp.tweak_plot_3D(fig,ax1,'general')
    rp.save_fig(outputfilename,fig,workdir)
