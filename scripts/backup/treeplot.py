#!/usr/bin/python2


#==============
def main():
#==============
    """
    Main function. Calls all the rest.
    """

    import treeplot_module as tm


    #-----------------------
    # Set up
    #-----------------------
    params = tm.global_params()
    params.read_cmdlineargs()
    params.get_output_info()

    if params.verbose:
        params.print_params()



    #---------------------------
    # read in data, make tree
    #---------------------------
    if params.from_backup:
        tree, outputnrs, t = tm.read_backup(params)

    else:
        descendants, progenitors, progenitor_outputnrs, outputnrs, t = tm.read_mergertree_data(params)
        tree = tm.make_tree(progenitors, descendants, progenitor_outputnrs, outputnrs, t, params)
        del progenitors
        del descendants
        del progenitor_outputnrs

        if params.write_backup:
            tm.dump_backup(tree, outputnrs, t, params)


    #------------------
    # Plot the tree
    #------------------
    tm.plot_tree(tree, outputnrs, t, params)

    #------------------------------
    # If needed: Plot particles
    #------------------------------
    if params.plotparticles:
        tm.plot_treeparticles(tree, t, params)


    return






#===============================
if __name__ == "__main__":
#===============================
    
    main()

