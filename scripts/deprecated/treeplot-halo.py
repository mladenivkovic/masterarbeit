#!/usr/bin/python3

#=======================================================
# Create merger tree plot for given halo at last output
#=======================================================



from os import getcwd
from sys import argv #command line arguments


import treeplot_module as tm
import ramsesplot_plot2D as rp




#===============================
if __name__ == "__main__":
#===============================

    #-----------------------
    # Set up
    #-----------------------

    workdir = str(getcwd())


    # Read args
    halo = int(argv[1])
    lastdir = str(argv[2])
    noutput = int(argv[3])
    ncpu= int(argv[4]) 

    # output filename for image
    outputfilename = "merger_tree_"+lastdir[-5:]+"_halo_"+str(halo)



    #-----------------------------
    # Print parameters to screen
    #-----------------------------

    print("===============================================")
    print("Working parameters are:")
    print("halo:", halo)
    print("last output directory:", lastdir)
    print("number of outputs:", noutput)
    print("ncpus used for sim:", ncpu)
    print("===============================================")


    #----------------
    # read in data
    #----------------
    descendants, progenitors, progenitor_outputnrs, outputnrs, t = tm.read_mergertree_data(noutput, ncpu, lastdir)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #--------------------
    # Debug mode:
    #--------------------

    # In case you want to test out, uncomment the following 8 lines. They contain
    # dummy arrays for testing. The data doesn't need to be read in anymore.

    #  progenitors = [ [202], [201,  150], [ 186, 263,   9], [182, 166], [   1,    2,   3,   4], [5,  6,  7, 8,  9, 10, 11, 12], [13, 14, 15, 16, 17, 18, 19, 20] ]
    #  descendants = [ [203], [202, -202], [-201, 201, 150], [186, 263], [-166, -182, 182, 166], [4, -2, -1, 1, -3,  2, -4,  3], [5,  6,  7, 8,  9, 10, 11, 12] ]
    #  progenitor_outputnrs = [ [7], [6, 6], [5, 5, 2], [4, 4], [3, 3, 3, 3], [2, 2, 2, 2, 2, 2, 2, 2], [1, 1, 1, 1, 1, 1, 1, 1]]
    #  noutput = len(progenitors)
    #  t = [ noutput-x for x in range(noutput) ]
    #  outputnrs = t
    #  halo= 203
    #  lastdir = 'output_00008'
    #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   

    # get number of last output
    lastdirnr = lastdir[-5:]
    has_leading_zeros = True
    while has_leading_zeros:
        has_leading_zeros = False
        if lastdirnr[0] == '0':
            lastdirnr = lastdirnr[1:]
            has_leading_zeros = True

    lastdirnr = int(lastdirnr)


    #------------------
    # Make tree
    #------------------
    tree = tm.make_tree(progenitors, descendants, progenitor_outputnrs, halo, lastdirnr)

    #------------------
    # Plot the tree
    #------------------
    fig, ax = tm.plot_tree(tree, outputnrs)

    #------------------
    # Tweak the plot
    #------------------
    tm.tweak_drawing(fig, ax, t, halo, lastdirnr)
    
    #------------------
    # Save the figure
    #------------------
    rp.save_fig(outputfilename, fig, workdir)


