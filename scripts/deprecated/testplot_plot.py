
# A module for the ramsesplot scripts
# Containts all functions about 2D plot

import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots_adjust
from matplotlib.font_manager import FontProperties # for legend



########################
### General variables
########################

explicit_colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', 
                    u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']




#====================================================
def plot_3D(x, y, z, ax, case, clumpID, colorindex):
#====================================================
    """
    plot given x,y,z. First gets parameters for case.
    x, y, z:    x, y, z coordinates to plot
    ax :        axis object
    case:       for which case to plot/choose parameters
    clumpID:    clump ID which is to be plotted.

    accepted cases:
        'child-mb'         most bound particles of child clump
        'child-transp'     child particles with low alpha 
        'halo-mb'          most bound particles of halo clump
        'halo-transp'      halo particles with low alpha 
        'level-halo'       halo particles for plotting by level
        'level-halo-mb'    most bound halo particles for plotting by level
        'level-subhalo'    subhalo particles for plotting by level
        'level-subhalo-mb' most bound subhalo particles for plotting by level


    returns:
        nothing

    """
    from mpl_toolkits.mplot3d import Axes3D

    mylabel, message, pointsize, pointalpha, mymarker, mylw, chose_colour, specific_colour  = get_plot_params(case, x.shape[0], clumpID, colorindex)

    if mylw == 0 :
        edgecolor = None
    else:
        edgecolor = u'#000000'

    print(message)

    if chose_colour:
        ax.scatter(x, y, z, 
                s=pointsize, 
                label=mylabel, 
                linewidth=mylw, 
                marker=mymarker, 
                depthshade=True, 
                alpha=pointalpha, 
                facecolor=specific_colour, 
                edgecolor=edgecolor)
    else:
        ax.scatter(x, y, z, 
                s=pointsize, 
                label=mylabel, 
                linewidths=mylw, 
                marker=mymarker, 
                depthshade=True, 
                alpha=pointalpha, 
                edgecolor=edgecolor)

    return







#===============================================================
def plot_32D(x, y, z, ax1, ax2, ax3, case, clumpID, colorindex):
#===============================================================
    """
    plot given x,y,z. First gets parameters for case.
    x, y, z:    x, y, z coordinates to plot
    ax :        axis object
    case:       for which case to plot/choose parameters
    clumpID:    clump ID which is to be plotted.
    colorindex: index to pick specific colour from

    accepted cases:
        'child-mb'         most bound particles of child clump
        'child-transp'     child particles with low alpha 
        'halo-mb'          most bound particles of halo clump
        'halo-transp'      halo particles with low alpha 
        'level-halo'       halo particles for plotting by level
        'level-halo-mb'    most bound halo particles for plotting by level
        'level-subhalo'    subhalo particles for plotting by level
        'level-subhalo-mb' most bound subhalo particles for plotting by level
        'seq-child'        child for sequence plotting
        'seq-halo'         halo for sequence plotting

    returns:
        nothing

    """

    mylabel, message, pointsize, pointalpha, mymarker, mylw, chose_colour, specific_colour  = get_plot_params(case, x.shape[0], clumpID, colorindex)
    print(message)


    if mylw == 0 :
        edgecolor = None
    else:
        edgecolor = u'#000000'


    if chose_colour:
        ax1.scatter(x, y, 
                s=pointsize, 
                label=mylabel, 
                marker=mymarker, 
                alpha=pointalpha, 
                edgecolor=edgecolor,
                linewidths=mylw, 
                c=specific_colour)
        ax2.scatter(y, z, 
                s=pointsize, 
                label=mylabel, 
                marker=mymarker, 
                alpha=pointalpha, 
                edgecolor=edgecolor,
                linewidths=mylw, 
                c=specific_colour)
        ax3.scatter(x, z, 
                s=pointsize, 
                label=mylabel, 
                marker=mymarker, 
                edgecolor=edgecolor,
                linewidths=mylw, 
                alpha=pointalpha, 
                c=specific_colour)
    else:
        ax1.scatter(x, y, 
                s=pointsize, 
                label=mylabel, 
                marker=mymarker, 
                edgecolor=edgecolor,
                linewidths=mylw, 
                alpha=pointalpha)
        ax2.scatter(y, z, 
                s=pointsize, 
                label=mylabel, 
                marker=mymarker, 
                edgecolor=edgecolor,
                linewidths=mylw, 
                alpha=pointalpha)
        ax3.scatter(x, z, 
                s=pointsize, 
                label=mylabel, 
                marker=mymarker, 
                edgecolor=edgecolor,
                linewidths=mylw, 
                alpha=pointalpha)


    return





#========================================================
def get_plot_params(case, length, clumpID, colorindex):
#========================================================
    """
    set up plot parameters for possible cases.
    case:      which case to plot.
               passed on by the plot_2D function.
    length:    lenght of array to plot
    clumpID:   clump ID
    colorindex: index to pick specific colour from

    returns:
        mylabel, message, pointsize, pointalpha, mymarker, mylw, chose_colour, specific_colour
    mylabel:        legend label for plot
    message:        message to print before plotting
    pointsize:      size of plot points
    pointalpha:     alpha of points to plot
    mymarker:       plot markers
    mylw:           line width around plot points
    chose_colour:   whether specific colour was specified
    specific_colour:specific colour chosen

    """



    #default values

    mylabel=None
    message='No message'
    pointsize=1
    mylw=0.0
    pointalpha=0.6
    mymarker='o'
    chose_colour = False
    specific_colour = ''




    if (case=='child-mb'):
        mylabel='Most bound of child ' + str(clumpID)
        message=mylabel +' : plotting ' + str(length) + ' particles'
        pointsize=2
        mylw=1
        pointalpha=0.8
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]
 

    elif (case=='child-transp'):
        mylabel='child ' + str(clumpID)
        message=mylabel +' : plotting ' + str(length) + ' particles'
        pointsize=2
        mylw=0.0
        pointalpha=0.2
        mymarker='.'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]


    elif (case=='halo-mb'):
        mylabel= 'Most bound of main '+str(clumpID)
        message= mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=2
        mylw=1.0
        pointalpha=0.6
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]



    elif (case=='halo-transp'):
        mylabel= 'Main '+str(clumpID)
        message= mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=2
        mylw=0.0
        pointalpha=0.1
        mymarker='.'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]



    elif (case=='level-halo'):
        mylabel= 'main subhalos '
        message= mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=1
        mylw=0.0
        pointalpha=0.2
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]



    elif (case=='level-halo-mb'):
        mylabel= 'main subhalos '
        message= mylabel + ' : plotting ' +str(length) + ' most bound particles'
        pointsize=2
        mylw=1.0
        pointalpha=0.6
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]


        
    elif (case=='level-subhalo'):
        mylabel= 'satellite subhalos level '+str(clumpID)
        message= mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=2
        mylw=0.0
        pointalpha=0.2
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]

        
    elif (case=='level-subhalo-mb'):
        mylabel= 'satellite subhalos level '+str(clumpID)
        message= mylabel + ' : plotting ' +str(length) + ' most bound particles'
        pointsize=2
        mylw=1.0
        pointalpha=0.6
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]



    elif (case=='seq-child'):
        mylabel= 'satellite halo '+str(clumpID)
        message= mylabel + ' : plotting ' +str(length) + '  particles'
        pointsize=4
        mylw=0.0
        pointalpha=0.8
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]


    elif (case=='seq-halo'):
        mylabel= 'main halo '+str(clumpID)
        message= mylabel + ' : plotting ' +str(length) + '  particles'
        pointsize=4
        mylw=0.0
        pointalpha=0.8
        mymarker='o'
        chose_colour = True
        specific_colour = explicit_colors[colorindex]


    else:
        print("get_plot_params: didn't recognize case. Quitting")
        quit()

   
    
    return mylabel, message, pointsize, pointalpha, mymarker, mylw, chose_colour, specific_colour 





