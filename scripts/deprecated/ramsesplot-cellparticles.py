#!/usr/bin/python3

# 
from os import getcwd
from sys import argv,exit #command line arguments

import ramsesplot_plot2D as rp
import ramsesplot_cells as rc 
import ramsesplot_fetch as rf




workdir= str(getcwd())

noutput=int(argv[1])
halo=int(argv[2])
nparticles=int(argv[3])

boxlen = 1000.0

outputfilename = "cell_particles_plot"+str(halo)



if __name__ == "__main__":


    # get particle data
    children, child_levels = rf.get_children_data(halo, noutput) 
    x_part, y_part, z_part, halox, haloy, haloz = rf.get_halo_particles(children, halo, noutput, nparticles)


    # get cell data
    cellpos, cell_info = rc.read_cell_data(4+2*noutput, noutput)
    x, y, z, levs = rc.filtercells(cellpos, cell_info[:,1], cell_info[:,0], halo)
    x_lev, y_lev, z_lev = rc.sortcells(x, y, z, levs)


    # plot cells
    fig, ax1, ax2, ax3 = rp.setup_3fig()
    rc.plotcells(x_lev, y_lev, z_lev, boxlen, fig, ax1, ax2, ax3)
    


    # plot particles
    rp.plot_32D(halox, haloy, haloz, ax1, ax2, ax3, 'cell-halo',halo)






    #tweak stuff
    rp.tweak_3plot_2D(fig, ax1, ax2, ax3, 'general')
    rp.save_fig(outputfilename, fig, workdir)


