#!/usr/bin/python3


#=======================================
# A module for handling mergertree data
#=======================================


#==================
class Node:
#==================
    """
    A class for each node to be drawn. 
    Parameters:
        id:         (sub)halo ID
        y:          y axis value for plot (time/redshift)
        desc_id:    the ID of the "main descendant". Needed to distinguish
                    jumps across multiple snapshots.
        
    """

    #-------------------------------------
    def __init__(self, ID, y, desc_id):
    #-------------------------------------
        """
        Arguments:
            id:         (sub)halo ID [any number]
            y:          y axis value for plot (time/redshift) [any number]
            desc_id:    the ID of the "main descendant". Needed to distinguish
                        jumps across multiple snapshots. [any number]
            
        """
        
        self.id = ID                # clump ID at timestep y
        self.main_desc_id = desc_id # id of "main descendant" to check for jumpers
        self.y = y                  # value for y axis: outputnumber where it was active clump
        self.x = 0                  # will be set later

        self.progs = []             # list of progenitor nodes
        self.is_main_prog = []      # wether the i-th progenitor of the progs list is the main progenitor

        self.branches = []          # list of all branches at this node
        self.branches_level = 0    # how many times this branch will split
        self.branches_tot = 0       # total number of branches for this branch
                                    # Needed for nice plotting
        self.start_of_branch = None # Node at which the branch this clump is in starts
        self.walked = False         # whether this node has been walked over already
        

        self.color = -1
        return



    #-----------------------------------------------
    def add_progenitor(self, prog_node, is_main):
    #-----------------------------------------------
        """
        Adds a progenitor to this node's list.
        Arguments:
            prog_node:  node of the progenitor to add to list [class Node object]
            is_main:    whether this progenitor is the main progenitor of this node [boolean] 
        """
        self.progs.append(prog_node)
        self.is_main_prog.append(is_main)
        return



#=========================================
def draw(node, ax, colorlist):
#=========================================
    """
    Plots all connections of Node node to its progenitors.
    First all main progenitors are plotted, going down straight
    lines. When the main progenitors are done, it starts plotting
    all the branches.

    Arguments:
        node:       class Node object whose progenitors are to be plotted
        ax:         axis object of the plot
        colorlist:  list of colors for different branches

    returns:
        nothing
    """


    
    # First go down straight lines along main progenitors
    for i, prog in enumerate(node.progs):
        if node.is_main_prog[i]:
            # call actual drawing function
            draw_on_plot(node, prog, ax, colorlist)
            # call yourself recursively for the main progenitor
            draw(prog, ax, colorlist)


    # Now go down branches
    for i, prog in enumerate(node.progs):
        if not node.is_main_prog[i]:
            # call actual drawing function
            draw_on_plot(node, prog, ax, colorlist)
            # call yourself recursively for this progenitor
            draw(prog, ax, colorlist)

    # if you reached the leaves, just add the progenitor numbers
    if len(node.progs) == 0:
        draw_on_plot(node, node, ax, colorlist)



    return

    



#====================================================
def draw_on_plot(node, prog, ax, colorlist):
#====================================================
    """
    The actual drawing function. Draws a line between Node node 
    and its progenitor prog.
    If the progenitor re-imerges at a later timestep, draw a dotted
    line instead of a solid line to signify where it looks like it 
    merged into.

    Arguments:
        node:       class Node object of descendant
        prog:       class Node object of progenitor
        ax:         axis object of the plot
        colorlist:  list of colors for different branches

    returns:
        nothing
    """

    # get x and y values for plot
    x = [node.x, prog.x]
    y = [node.y, prog.y]


    # get line colors
    outer = prog.color//len(colorlist)
    inner = prog.color - outer * len(colorlist)
    outer = -(outer+1)

    # Determine line-style
    # Dashed for mergers that will re-emerge later
    linestyle='-'
    linewidth = 6
    outerlinewidth = 10
    alp = 1

    if node.id != prog.id: 
        if node.id != prog.main_desc_id: 
            linestyle = '--'
            linewidth = 3
            outerlinewidth = 4
            alp = 0.2



    #---------------
    # Plot the line
    #---------------

    # plot outer line
    ax.plot(x, y,
            color = colorlist[outer],
            zorder=1,
            lw=outerlinewidth,
            ls='-',
            alpha=alp)

    # plot inner line
    ax.plot(x, y, 
            color = colorlist[inner], 
            zorder=2, 
            lw=linewidth,
            ls=linestyle)


    #------------------------
    # Plot points at ends
    #------------------------

    #  ax.scatter(node.x, node.y,
    #          c=colorlist[inner],
    #          linewidth=3,
    #          edgecolor=colorlist[outer],
    #          marker='o',
    #          s=150,
    #          zorder=3)



    #---------------------
    # Annotation
    #---------------------


    # plot line from scatterpoint to annotation
    #  x = [node.x, node.x+xoffset]
    #  y = [node.y, node.y+yoffset]

    # outer line
    #  linestyle = '-'
    #  linewidth = 2
    #  ax.plot(x, y,
    #          ls = linestyle,
    #          lw = linewidth,
    #          color = colorlist[outer],
    #          zorder = 1)


    # inner line
    #  linestyle = ':'
    #  linewidth = 2
    #  ax.plot(x, y,
    #          ls = linestyle,
    #          lw = linewidth,
    #          color = colorlist[inner],
    #          zorder = 1)

    # take here node data instead of prog!
    outer = node.color//len(colorlist)
    inner = node.color - outer * len(colorlist)
    outer = -(outer+1)


    # Annotate the dots with the clump ID

    bbox_props = dict(boxstyle="round,pad=0.1", 
            fc=colorlist[inner], 
            ec=colorlist[outer], 
            lw=2, 
            alpha=1)

    t = ax.text(node.x, node.y, str(node.id),
            size=10,
            bbox=bbox_props,
            horizontalalignment = 'center',
            verticalalignment = 'center')




    return





#=====================================
def get_x(node, colorindex, borders):
#=====================================
    """
    Assign x values and colors for the plot to a node and its progenitors.
    First descend down branches to determine the space required and
    mark the nodes of the branch as walked over, so you get nice
    straight lines for multisnapshot-jumpers.
    For main progenitors, just inherit the x and color values.



    Arguments:
        node:       class Node object to check for
        colorindex: currently last used colorindex
        borders:    list for [xmin, xmax] used in tree. Needed for plotting later.

    returns:
        colorindex: currently last used colorindex
        borders:    current borders

    """

    # First descend down branches, not straight lines,
    # to fix the x coordinates of multisnapshot jumps
    for i, prog in enumerate(node.progs):
        if (not node.is_main_prog[i]) and (not prog.walked):

            # if there is more than one merger at this y, add 
            # previously determined branchindex

            ids = (branch.id for branch in prog.start_of_branch.branches)
            ids = list(ids)
            ind = ids.index(prog.id)


            # figure out whether to go left or right
            if  prog.start_of_branch.x > prog.start_of_branch.start_of_branch.x:
                # go right first
                sig = 1
            else:
                # go left first
                sig = -1

            if (len(ids)-ind)%2 == 0:
                sig = -sig


            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # hic sunt dracones
            skip = 1
            for b in range(ind+1):
                branch = prog.start_of_branch.branches[b]
                if prog.id != branch.id:
                    skip += (branch.branches_level+1) * branch.branches_tot + 1
                else:
                    skip += (branch.branches_level+1) * branch.branches_tot/2 + 1
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


            # set new value for new branch
            prog.x = node.x + sig * skip

            # update borders
            if prog.x > borders[1]:
                borders[1] = prog.x
            if prog.x < borders[0]:
                borders[0] = prog.x

            # give new branch a new color
            colorindex += 1
            prog.color = colorindex
            prog.walked = True


            # call yourself recursively
            colorindex, borders = get_x(prog, colorindex, borders)


 
    # if it's its main prog, inherit x and color
    for i, prog in enumerate(node.progs):
        if node.is_main_prog[i] and (not prog.walked):
            prog.x = node.x
            prog.color = node.color
            prog.walked = True
            colorindex, borders = get_x(prog, colorindex, borders)




    return colorindex, borders





#==============================================================================
def make_tree(progenitors, descendants, progenitor_outputnrs, halo, lastdirnr):
#==============================================================================
    """
    makes a tree out of read in lists. 
    Thre tree is a list containing lists of Nodes (see class Node above)
    for each output step.

    parameters:
        progenitors:            list of lists of progenitors for each timestep
        descendants:            list of lists of descendants for each timestep
        progenitor_outputnrs:   list of lists of the output numbers when 
                                any given progenitor was an active clump
        halo:                   root halo ID
        lastdirnr:              the number of the last output_XXXXX directory

    returns:
        tree:                   list of lists of nodes constituting the tree
    """

    
    #------------------
    # Setup
    #------------------

    nout_present = len(progenitors) 
    output_start = lastdirnr - nout_present



    #---------------------
    # initialise tree
    #---------------------

    tree = []

    # create empty list for each output
    for i in range(nout_present+1):
        tree.append([])

    # enter root
    rootnode = Node(halo, lastdirnr, 0)
    tree[0]=[rootnode]

    

    #---------------------
    # Make tree
    #---------------------

    print()

    # Loop over all snapshots
    for out in range(nout_present):
        # for each branch of the tree at that snapshots:
        for branch in tree[out]:

            # find for which progenitor the descendant matches
            for i in range(len(progenitors[out])) :

                snapnr = progenitor_outputnrs[out][i]       # snapshot nr for progenitor
                ind = nout_present + output_start - snapnr  # index in tree / tree level
                progid = abs(progenitors[out][i])           # progenitor ID

                if abs(descendants[out][i]) == branch.id:

                    is_main = descendants[out][i] == branch.id
                    # is a main progenitor if desc ID == branch ID
                    # is a merger if desc ID == -branch ID

                    # Check first if progenitor is already in list
                    prog_not_in_list = True

                    if (progid > 0):    # always add a new 0!
                        if (len(tree[ind]) > 0):
                            for j, candidate in enumerate(tree[ind]):
                                if (candidate.id == progid):
                                    # Then this progenitor is already in the tree.
                                    branch.add_progenitor(tree[ind][j], is_main)
                                    prog_not_in_list = False
                                    break

                    if prog_not_in_list:
                        # create new node for progenitor
                        newnode = Node(progid, snapnr, branch.id)

                        # add new node to tree
                        tree[ind].append(newnode)

                        # link to its progenitor:
                        # you know in which outputnr it will be
                        # since it's the last added element, it's index will
                        # be len(tree at that outputnr) - 1
                        branch.add_progenitor(tree[ind][len(tree[ind])-1], is_main)


                    # Print informations
                    print('Adding progenitor ', end=' ')
                    print('{0:6d}'.format(progid), end=' ')
                    print('for descendant ', end=' ' )
                    print('{0:6d}'.format(branch.id), end=' ')
                    print("|| snapshot " , end=' ' )
                    print('{0:3d}'.format(snapnr), end=' ')
                    print("->" , end=' ' )
                    print('{0:3d}'.format(branch.y), end=' ')



                    if (is_main):
                        if branch.y-snapnr > 1 :
                            print("   detected jumper")
                        else:
                            print() # make newline
                    else:
                        if branch.y-snapnr > 1 :
                            print("   detected jumper and merger")
                        else:
                            print("   detected merger")


        print() # make new line after every snapshot read-in


    # remove empty lists at the end of the tree:
    # clumps might not have been around since the first output

    list_is_empty = (len(tree[-1]) == 0)

    while list_is_empty:
        del tree[-1]
        list_is_empty = len(tree[-1]) == 0

            


    return tree




#=============================
def plot_tree(tree, yaxis):
#=============================
    """
    The main function for plotting the tree.

    Arguments:
        tree:   list of lists of class Node objects, constituting the tree
        yaxis:  labels for the (left) y-Axis

    returns:
        fig:    pyplot figure object
        ax:     axis object of plot
    """

    import matplotlib.pyplot as plt

    # Define a long list of colors for the branches
    longcolorlist=[ 'red', 'green', 'gold', 'magenta', 'orange', 'mediumpurple', 'cyan', 'lime', 'saddlebrown', 'lightpink','mediumseagreen']



    # First find the number of branches for each branch.
    tree[0][0].start_of_branch = tree[0][0]
    walk_tree(tree[0][0], tree[0][0])

    # Now recursively sum up the number of branches
    # to know how much space to leave between them for plotting
    sum_branches(tree[0][0])
    sort_branch(tree[0][0])


    # reset whether nodes have been walked
    for level in tree:
        for branch in level:
            branch.walked = False


    # start distributing x values and colors. Set root at x = 0
    tree[0][0].x = 0
    colorindex = 0
    tree[0][0].color = colorindex


    borders = [0, 0, tree[0][0].y, tree[-1][0].y]

    ignore_this_output, borders = get_x(tree[0][0], colorindex, borders)

    dx = borders[1]-borders[0]
    if (dx == 0):
        dx = 1
    borders = [borders[0]-dx/20, borders[1]+dx/20]



    print()
    print("Creating figure.")
    print()

    # create figure
    fig = plt.figure()
    fig.set_size_inches(40, 20)
    ax = fig.add_subplot(111)

    ax.set_xlim(borders)

    # draw the tree
    draw(tree[0][0], ax, longcolorlist)



    return fig, ax





#==================================================
def read_mergertree_data(noutput, ncpu, lastdir):
#==================================================
    """
    reads in mergertree data as written by the mergertree patch.
    parameters:
        noutput:    number of output directories to go through
        ncpu:       number of processors to go through
        lastdir:    last directory of simulation run, i.e. 
                    directory to start mergertree with

    returns:
        progenitors :           lists of lists of progenitors and descendants,
        descendants :           starting with the last output step.
        progenitor_outputnrs:   the output number at which the progenitor is
        outputnrs   :           the output number at which descendants were taken from
        time :                  list of times correspondig to each output step

    """ 

    import numpy as np
    import warnings


    print("Reading in mergertree data.")

    


    # create lists where to store stuff
    fname = 'mergertree.txt'

    progenitors = []
    descendants = []
    time = []
    progenitor_outputnrs = []
    outputnrs = []

    startnr=int(lastdir[-5:])
    dir_template = lastdir[:-5] # = 'output_'


    # Loop over directories
    for output in range(noutput - 1): #-1: No output from first directory needed!!!!
        # Start with last directory (e.g. output_00060),
        # work your way to first directory (e.g. output_00001)
        dirnr =  str(startnr - output).zfill(5)
        srcdir = dir_template + dirnr
        progs_snapshot = []
        descs_snapshot = []
        prog_outputnr_snapshot = []

        # Stop early if you reach a directory that has no mergertree.txt* files
        # (Can happen if there are no halos in the simulation yet)
        try:
            # loop over files
            for cpu in range(ncpu):
                filenr = str(cpu + 1).zfill(5)          # get 00001, 00002, ...
                fileloc = srcdir + '/' + fname + filenr


                # ignore "empty files" warnings
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    raw_data = np.loadtxt(fileloc, dtype='int', skiprows=1, usecols=([0, 1, 2]))

                #intercept index error if there is only one halo
                try:
                    dshape = raw_data.shape[1] #try if the shape is 2dimensional

                    desc = raw_data[:,0].tolist()
                    prog = raw_data[:,1].tolist()
                    prog_outputnr = raw_data[:,2].tolist()

                except IndexError: 

                    # in this case, the read file has only 1 halo, or none.
                    # Check if there is 1:

                    if (len(raw_data)>0):
                        desc = [raw_data[0]]
                        prog = [raw_data[1]]
                        prog_outputnr = [raw_data[2]]
                    else:
                        continue

                # append read list for each CPU to list for entire snapshot
                for i in range(len(prog)):
                    progs_snapshot.append(prog[i])
                    descs_snapshot.append(desc[i])
                    prog_outputnr_snapshot.append(prog_outputnr[i])
            



            #  print("==============================================================")
            #  print("READCHECK:", dirnr, " desc", descs_snapshot, "prog", progs_snapshot)
            #  print("==============================================================")
            #  print()

            # append entire list here!
            descendants.append(descs_snapshot)
            progenitors.append(progs_snapshot)
            progenitor_outputnrs.append(prog_outputnr_snapshot)
            outputnrs.append(startnr - output)


            # get time
            fileloc = srcdir+'/info_'+dirnr+'.txt'
            infofile = open(fileloc)
            for i in range(8):
                infofile.readline() # skip first 8 lines
            
            timeline = infofile.readline()
            timestring, equal, timeval = timeline.partition("=")
            timefloat = float(timeval)

            time.append(timefloat)
    
        except OSError: # If file doesn't exist
            break


    return descendants, progenitors, progenitor_outputnrs, outputnrs, time
   




#===========================
def sort_branch(branch):
#===========================
    """
    Sort the list of branches of a given start of a branch (class Node object)
    by increasing time of appearance. The earliest branches (lowest on y-axis)
    come first. If there are multiple mergers at the same time, put the one 
    with more own branches further out.

    Parameters:
        branch:     class Node objecs whose .branches list is to be sorted

    returns:
        nothing
    """
            
    from copy import deepcopy

    if len(branch.branches) > 0:
        branchlist = (bbranch for bbranch in branch.branches)
        branchlist = list(branchlist)
        branchlist_copy = deepcopy(branchlist)

        times = []
        all_branches = []

        for bbranch in branchlist:
            times.append(bbranch.y)
            all_branches.append(bbranch.branches_tot)

        needs_branchcount_sorting = False

        for t in times:
            if times.count(t) > 1:
                needs_branchcount_sorting = True
                break

        sort_ind = range(len(times))
        times, sort_ind = (list(i) for i in zip(*sorted(zip(times, sort_ind))))

        for i in range(len(times)):
            branchlist[i] = branchlist_copy[sort_ind[i]]

        # make a copy of sorted list
        branchlist_copy = deepcopy(branchlist)

        if needs_branchcount_sorting:
            # times[] is now sorted. check whether
            # at least the next element is the same
            i = 0
            while i < (len(times)-1):
                if (times[i+1]==times[i]):
                    # if the following is the same:
                    startind = i
                    endind = i+1

                    # find where the same end
                    while times[startind]==times[endind]:
                        endind += 1
                        if endind == len(times):
                            break

                    # sort all of them by ascending branch_tot
                    branchtots = []
                    sort_ind2 = []
                    for j in range(startind, endind-1):
                        branchtots.append(all_branches[sort_ind[j]])
                        sort_ind2.append(j)

                    branchtots, sort_ind2 = (list(b) for b in zip(*sorted(zip(branchtots, sort_ind2))))

                    for k in range(endind-startind-1):
                        branchlist[startind+k] = branchlist_copy[sort_ind2[k]]


                    i = endind

                else:
                    i+=1


        #overwrite branche's branch list
        for i in range(len(times)):
            branch.branches[i] = branchlist[i]




    return





#=================================
def sum_branches(node):
#=================================
    """
    Recursively sum up the total number of branches of each "root"
    to effectively determine the space needed for decent plotting.

    Arguments:
        node:   class Node object to check for

    returns:
        nothing
    """


    # each branch root has itself as a branch.
    # If the tree splits up somewhere else along the way,
    # then there must be > 1 branch in the root node.
    
    if len(node.branches) > 0:
        node.branches_level = 1
        for i, branch in enumerate(node.branches):
            # don't call yourself again
            sum_branches(branch)
            node.branches_tot += branch.branches_tot
            node.branches_tot += 1 # count yourself too :)
            
            node.branches_level += branch.branches_level
            sort_branch(branch)






    return





#==================================================
def tweak_drawing(fig, ax, yaxis, halo, lastdirnr):
#==================================================
    """
    tweaks the plot. Removes x-axis ticks, labels y-axis
    ticks nicely, adds right y axis.

    parameters:
        fig:        pyplot figure object
        ax:         axis object of plot
        yaxis:      array for y axis ticks
        halo:       halo number needed for title
        lastdirnr:  number of last directory, where the root of the tree starts

    returns:
        nothing
    """


    import matplotlib.pyplot as plt

    #-----------------
    # preparation
    #-----------------

    yaxis = [yaxis[0]+1] + yaxis # add halo output step to yaxis
    noutput = len(yaxis)
    firstoutput = lastdirnr - noutput + 2
    outputnr = range(lastdirnr, firstoutput-2, -1)
    # firstoutput -2:
    # -1 because range makes list [from, to)
    # -1 because you're plotting progenitor at each step, which is at a timestep earlier







    #---------------------
    #  Prepare left y axis
    #---------------------

    # determine how many y axis ticks you want.
    # find step to go through loop

    nyticks_step = int(noutput/10) + 1
    yticks = []
    
    ind = 0
    while ind < noutput:
        
        if ind % nyticks_step == 0:
            yticks.append(outputnr[ind])
    
        ind += 1
    
    


    #----------------------
    # Prepare right y axis
    #----------------------

    ax2 = ax.twinx()
    ax2.set_ylim([firstoutput-2,lastdirnr+1])
    
    
    yticks_right = []
    yticks_right_labels=[]
    
    ind = 0
    while ind < noutput:
        if ind % nyticks_step == 0:
            yticks_right.append(outputnr[ind])
            yticks_right_labels.append(round(yaxis[ind],2))

        ind += 1




    #---------------------
    # Set ticks and title
    #---------------------

    # left y ticks
    ax.set_yticks(yticks)
    ax.set_ylim([firstoutput-2, lastdirnr+1])
    ax.set_ylabel('output number', size=20)
    ax.tick_params(axis='both', labelsize=15)

    # right y ticks
    ax2.set_yticks(yticks_right)
    ax2.set_yticklabels(yticks_right_labels)
    ax2.set_ylabel("t [code units]", size=20)
    ax2.tick_params(axis='both', labelsize=15)


    # x ticks
    ax.set_xticks([]) # no x ticks

    # title
    title = "Merger tree for halo "+str(halo)+" at output "+str(lastdirnr)
    ax.set_title(title, size=26)

    



    # add grid
    ax.grid()

    # cleaner layout
    plt.tight_layout()
    
    return





#=================================
def walk_tree(node, root):
#=================================
    """ 
    Walk the tree and count the branches. Add branches to the 
    root/start of that branch.

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
                root.branches.append(prog)
                prog.start_of_branch = root

                walk_tree(prog, prog)


    # then just resume where you left off
    for i, prog in enumerate(node.progs):
        if not prog.walked:
            if node.is_main_prog[i]:
                prog.start_of_branch = root
                walk_tree(prog, root)

    return


