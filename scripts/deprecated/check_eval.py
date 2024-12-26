#!/usr/bin/python2


#========================================================================
# this script computes the halo properties at z=0
# for my local fieldtest.
# It does it in a way more inefficient way, but
# i know it to work.
#
#
#
# usage:
# check_eval.py <some integer > 0 >
#
#
# it imports a few functions from treeplot_module.py, which
# I didnt want to change. That's why it needs an integer
# (in treeplot, that is the required halo number to work with)
#
# Note: Like treeplot.py, this script finds the directory with z ~ 0
# by itself. If you want another directory, use -s <dirnr>
#========================================================================


import treeplot_module as tm
import numpy as np




#===================================
def eval_tree(tree, halo, params):
#===================================

    """
    Main function for tree evaluation.
    Parameters:
        tree:   tree as created by make_tree(...)
        halo:   halo ID for which to work
        params: global params object

    returns: 
        branchlen:  length of main branch
        nbranches:  number of branches in tree (excluding main branch)
    """

    # First find the number of branches for each branch.
    walk_tree(tree[0][0], tree[0][0])
    branchlen = find_main_branch_length(tree[0][0], 0)
    
    nbranches=tree[0][0].branches_tot

    return branchlen, nbranches




#=================================================
def find_main_branch_length(node, branchlength):
#=================================================
    """ 
    Walk down the main branch only using recursion, get branchlength.

    Parameters:
        node:           current node in the tree
        branchlength:   current branchlength

    returns:
        branchlength:   updated branchlength. 
    """

    # First check out new branches to mark possible jumps
    # over multiple timesteps as "walked"

    for i, prog in enumerate(node.progs):
        if node.is_main_prog[i]:
            branchlength = find_main_branch_length(prog, branchlength)
            branchlength += node.y - prog.y
            break


    return branchlength






#=================================
def walk_tree(node, root):
#=================================
    """ 
    Walk the tree and count the branches. Add branches to the 
    root/start of that branch in root.branches_tot

    Arguments:
        node:   class Node object to check whether a new branch starts here
        root:   class Node object which is current source of the branch.

    returns:
        nothing
    """

    # mark node as walked.
    node.walked = True

    # First check out new branches to mark possible jumps
    # over multiple timesteps as "walked"

    for i, prog in enumerate(node.progs):
        # check only if node hasn't been walked over already:
        if not prog.walked:
            if not node.is_main_prog[i]:
                # a new branch starts!
                # this progenitor will be the root for the new branch.
                root.branches_tot += 1
                walk_tree(prog, root)


    # then just resume where you left off
    for i, prog in enumerate(node.progs):
        if not prog.walked:
            if node.is_main_prog[i]:
                walk_tree(prog, root)

    return






#======================
def main():
#======================

    #------------------------
    # Read in stuff
    #------------------------

    params = tm.global_params()
    params.read_cmdlineargs()
    params.get_output_info()
    params.set_outputfilename()

    if params.verbose:
        params.print_params()

    descendants, progenitors, progenitor_outputnrs, outputnrs, t = tm.read_mergertree_data(params)

    #------------------------------------
    # Find which output to start with
    #------------------------------------
    startind = 0
    if params.start > 0:
        startind = np.where(outputnrs == params.start)[0][0]
    else:
        if not params.use_t:
            # find output closest to z=0
            startind = np.argmin(np.absolute(t))

    print "starting output is", outputnrs[startind]

    #-------------------------------------------------------------
    # Loop through all haloes in starting snapshot and evaluate
    #-------------------------------------------------------------
    for i,halo in enumerate(descendants[startind]):
        if halo > 0:

            # reset global params for new tree
            params.halo = halo
            params.start = 0

            tree = tm.make_tree(progenitors, descendants, progenitor_outputnrs, outputnrs, t, params)

            branchlen, nbranches = eval_tree(tree, halo, params)
            del tree

            print i, "halo", halo, "has length of main branch",
            print branchlen, " and ", nbranches, "branches"



    return






#===================================
if __name__ == "__main__":
#===================================
    
    main()
