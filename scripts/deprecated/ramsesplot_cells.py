#!/usr/bin/python3


#====================================================================
# This module contains functions needed for the plotting/working with
# Cells.
#====================================================================

# 
from sys import argv,exit #command line arguments

import warnings

import matplotlib 
# matplotlib.use('Agg') #don't show anything unless I ask you to. So no need to get graphical all over ssh.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties # for legend
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size



def read_cell_data(start_ind, noutput):
    """
    reads in cell data from files.
    start_ind:      index from which to start out from when reading filenames 
                    given as cmd line args
    noutput:        number of files to go through
    
    returns:
        data:     numpy array of the data
            data[0]: x coordinate of cell
            data[1]: y coordinate of cell
            data[2]: z coordinate of cell
        data_int: numpy array of integer data
            data_int[0]: cell level
            data_int[1]: cell clump ID
    """


    print("Reading in cell data.")


    # Loop over files
    for i in range(0,noutput):
        inputfile=str(argv[i+start_ind]) 
        
        # get cell center
        # ignore empty files
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            temp_data=np.loadtxt(inputfile, dtype='float', skiprows=1, usecols=[0,1,2], delimiter=',')
        if(temp_data.shape[0]>0):
            if 'data' in locals():
                data = np.vstack((data, temp_data))
            else:
                data = temp_data
        
        #get clump ids, parents and levels
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            temp_data_int=np.loadtxt(inputfile, dtype='int', skiprows=1, usecols=[3,4], delimiter=',')
        if(temp_data_int.shape[0]>0):
            if 'data_int' in locals():
                data_int= np.vstack((data_int, temp_data_int))
            else:
                data_int = temp_data_int


    return  data, data_int





def filtercells(cellpos, IDs,levels,clumpID):
    """
    filters cells that belong to given ID.
    cellpos:    3-dim array of cell positions
    IDs:        array of peak patch IDs of cells
    levels:     array of cell levels
    clumpID:    clumpID for which to filter for


    returns:
        x, y, z :   list of filtered cell positions
        levs:       list of filtered cell levels
    """

    ncells = 0
    for cell in IDs:
        if cell == clumpID:
            ncells += 1

    # x = np.zeros(ncells)
    # y = np.zeros(ncells)
    # z = np.zeros(ncells)
    # levs = np.zeros(ncells)
    x = [0]*ncells
    y = [0]*ncells
    z = [0]*ncells
    levs = [0]*ncells

    position = 0
    for i in range(len(IDs)):
        if IDs[i] == clumpID:
            x[position]=cellpos[i][0]
            y[position]=cellpos[i][1]
            z[position]=cellpos[i][2]
            levs[position]=levels[i]
            position += 1


    return x, y, z, levs




def sortcells(x,y,z,levs):
    """
    Sort cells by level.
    x: x position of cells
    y: y position of cells
    z: z position of cells

    returns:
        x_lev, y_lev, z_lev :   list of lists of cells by level.
                                x_lev[3] are cells of level 3.
    """
    maxlevel = int(max(levs))
    print ("maxlevel is:", maxlevel)
    ncells = len(x)

    cellcounter = [0]*(maxlevel+1)

    #count how many cells there are for each level
    for cell in levs:
        cellcounter[cell] += 1

    x_lev = []
    y_lev = []
    z_lev = []
    for i in range(maxlevel+1):
        x_lev.append(np.zeros(cellcounter[i])) 
        y_lev.append(np.zeros(cellcounter[i])) 
        z_lev.append(np.zeros(cellcounter[i])) 

    #reset cell counter
    cellcounter = [0]*(maxlevel + 1)

    #fill up arrays
    for i, cell in enumerate(levs):
        x_lev[cell][cellcounter[cell]] = x[i]
        y_lev[cell][cellcounter[cell]] = y[i]
        z_lev[cell][cellcounter[cell]] = z[i]
        cellcounter[cell] += 1


    return x_lev, y_lev, z_lev 




def plotcells(x_lev, y_lev, z_lev, boxlen, fig, ax1, ax2, ax3):
    """
    Plots cells according to their level.
    x_lev, y_lev, z_lev: list of numpy arrays of cells by level
    boxlen:              boxlen
    fig:                 figure object
    ax1, ax2, ax3:       axis object to draw on
    """

    maxlevel = len(x_lev)

    plotted = False

    cellalpha = 0.5

    for lev in range(maxlevel):
        if len(x_lev[lev]) > 0: #if there are cells of this level
            print("Plotting", x_lev[lev].shape[0], "cells of level", lev )
            scat1 = ax1.scatter(x_lev[lev], y_lev[lev], label="cells of level "+str(lev), marker='s', alpha = cellalpha)
            scat2 = ax2.scatter(y_lev[lev], z_lev[lev], label="cells of level "+str(lev), marker='s', alpha = cellalpha)
            scat3 = ax3.scatter(x_lev[lev], z_lev[lev], label="cells of level "+str(lev), marker='s', alpha = cellalpha)
            
            fig.canvas.draw()
            # if not plotted:
            #     fig.canvas.draw()
            #     plotted = True

        #     # Calculate cell size in pixels :
            size = boxlen/2**lev
            cs_pix1 = (ax1.transData.transform(np.vstack([size, size]).T) - ax1.transData.transform(np.vstack([np.zeros(1), np.zeros(1)]).T))
            cs_pix2 = (ax2.transData.transform(np.vstack([size, size]).T) - ax2.transData.transform(np.vstack([np.zeros(1), np.zeros(1)]).T))
            cs_pix3 = (ax3.transData.transform(np.vstack([size, size]).T) - ax3.transData.transform(np.vstack([np.zeros(1), np.zeros(1)]).T))
            spix1, _ = cs_pix1.T
            spix2, _ = cs_pix2.T
            spix3, _ = cs_pix3.T
        # 
         # 
        #     # Calculate and update size in points:
            size_pt1 = (spix1/fig.dpi*72)**2
            size_pt2 = (spix2/fig.dpi*72)**2
            size_pt3 = (spix3/fig.dpi*72)**2
            scat1.set_sizes(size_pt1)
            scat2.set_sizes(size_pt2)
            scat3.set_sizes(size_pt3)


    return





