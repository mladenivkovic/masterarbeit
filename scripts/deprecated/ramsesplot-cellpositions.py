#!/usr/bin/python3

# 
from os import getcwd
from sys import argv,exit #command line arguments

import ramsesplot_plot2D as rp
import ramsesplot_cells as rc 




workdir= str(getcwd())

noutput=int(argv[1])
halo=int(argv[2])
particles=int(argv[3])

boxlen = 1000.0

outputfilename = "cell_positions_plot"+str(halo)



if __name__ == "__main__":

    # get data

    cellpos, cell_info = rc.read_cell_data(4, noutput)
    x, y, z, levs = rc.filtercells(cellpos, cell_info[:,1], cell_info[:,0], halo)
    x_lev, y_lev, z_lev = rc.sortcells(x, y, z, levs)

    fig, ax1, ax2, ax3 = rp.setup_3fig()

    rc.plotcells(x_lev, y_lev, z_lev, boxlen, fig, ax1, ax2, ax3)
    
    rp.tweak_3plot_2D(fig, ax1, ax2, ax3, 'general')
    rp.save_fig(outputfilename, fig, workdir)


