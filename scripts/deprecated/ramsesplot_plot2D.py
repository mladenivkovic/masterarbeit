
# A module for the ramsesplot scripts
# Containts all functions about 2D plot

import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots_adjust
from matplotlib.font_manager import FontProperties # for legend





#====================
def setup_fig():
#====================
    """
    creates a new figure object.
    returns:
        fig, ax

        fig: figure object
        ax:  axis object (subplot)


    """
    print("Creating figure")
    fig = plt.figure(facecolor='white', figsize=(10,11))
    ax = fig.add_subplot(111)


    return fig, ax 



#====================
def setup_3fig():
#====================
    """
    creates a new figure object with 3 subplots.
    returns:
        fig, ax1, ax2, ax3

        fig: figure object
        ax123:  axis object (subplot)


    """
    print("Creating figure")
    fig = plt.figure(facecolor='white', figsize=(30,11))
    ax1 = fig.add_subplot(131, aspect='equal')
    ax2 = fig.add_subplot(132, aspect='equal')
    ax3 = fig.add_subplot(133, aspect='equal')


    return fig, ax1, ax2, ax3





#============================================
def save_fig(this_name,fig,workdir):
#============================================
    """
    Save figure as png.
    this_name:  name to save figure with
    fig:        figure object
    workdir:    string of directory where to save figure

    returns:
        nothing

    """

    fig_path = workdir+'/'+this_name+'.png'
    print("saving figure as "+fig_path)
    plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=100)#,bbox_inches='tight' )
    print("\n Saved", this_name+'.png')


    plt.close()
   
    return




#============================================
def plot_2D(x, y, z, ax, case, clumpID):
#============================================
    """
    plot given x,y,z. First gets parameters for case.
    x, y, z:    x, y, z coordinates to plot
    ax :        axis object
    case:       for which case to plot/choose parameters
    clumpID:    clump ID which is to be plotted.


    accepted cases:
        'all'           plot all particles that are in a halo
        'cell-halo'     plot all particles in halo for cellpart plot
        'child'         child clump particles
        'halo'          halo namegiver particles
        'level-halo'    halo for plotting clumps by level
        'level-subhalo' subhalo for plotting clumps by level

    returns:
        nothing

    """

    mylabel, message, pointsize, pointalpha, mymarker, mylw, chose_colour, specific_colour = get_plot_params(case, x.shape[0], clumpID)

    print(message)

    if chose_colour:
        ax.scatter(x, y, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw, 
                marker=mymarker, 
                alpha=pointalpha, 
                c=specific_colour)
    else:
        ax.scatter(x, y, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw, 
                marker=mymarker, 
                alpha=pointalpha)

    return



#==================================================
def plot_32D(x, y, z, ax1, ax2, ax3, case, clumpID):
#==================================================
    """
    plot given x,y,z. First gets parameters for case.
    x, y, z:        x, y, z coordinates to plot
    ax1,ax2,ax3:    axis object
    case:           for which case to plot/choose parameters
    clumpID:        clump ID which is to be plotted.


    accepted cases:
        'all'           plot all particles that are in a halo
        'cell-halo'     plot all particles in halo for cellpart plot
        'child'         child clump particles
        'halo'          halo namegiver particles
        'level-halo'    halo for plotting clumps by level
        'level-subhalo' subhalo for plotting clumps by level

    returns:
        nothing

    """


    mylabel, message, pointsize, pointalpha, mymarker, mylw, chose_colour, specific_colour = get_plot_params(case, x.shape[0], clumpID)

    print(message)

    if chose_colour:
        ax1.scatter(x, y, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw,
                marker=mymarker, 
                alpha=pointalpha, 
                c=specific_colour)
        ax2.scatter(y, z, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw, 
                marker=mymarker, 
                alpha=pointalpha, 
                c=specific_colour)
        ax3.scatter(x, z, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw, 
                marker=mymarker, 
                alpha=pointalpha, 
                c=specific_colour)
    else:
        ax1.scatter(x, y, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw, 
                marker=mymarker, 
                alpha=pointalpha)
        ax2.scatter(y, z, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw, 
                marker=mymarker, 
                alpha=pointalpha)
        ax3.scatter(x, z, 
                s=pointsize, 
                label=mylabel, 
                lw=mylw, 
                marker=mymarker, 
                alpha=pointalpha)


    return





#============================================
def get_plot_params(case, length, clumpID):
#============================================
    """
    set up plot parameters for possible cases.
    case:      which case to plot.
               passed on by the plot_2D function.
    length:    lenght of array to plot
    clumpID:   clump ID

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


    if (case=='all'):
        mylabel='halo particles'
        message=mylabel +' : plotting ' + str(length) + ' particles'
        pointsize=1
        mylw=0.0
        pointalpha=1
        mymarker=','


    elif (case=='cell-halo'):
        mylabel= 'halo '+str(clumpID) + " particles"
        message= mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=2
        mylw=1
        pointalpha=.5
        mymarker='o'
        chose_colour = True
        specific_colour = 'r'


    elif (case=='child'):
        mylabel='child ' + str(clumpID)
        message=mylabel +' : plotting ' + str(length) + ' particles'
        pointsize=2
        mylw=0.0
        pointalpha=0.6
        mymarker='o'


    elif (case=='halo'):
        mylabel= 'namegiver '+str(clumpID)
        message='halo ' + mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=1
        mylw=0.0
        pointalpha=0.6
        mymarker='o'



    elif (case=='level-halo'):
        mylabel= 'main subhalos '
        message= mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=1
        mylw=0.0
        pointalpha=0.6
        mymarker='o'

        
    elif (case=='level-subhalo'):
        mylabel= 'satellite subhalos level '+str(clumpID)
        message= mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=1
        mylw=0.0
        pointalpha=0.6
        mymarker='o'

    else:
        print("get_plot_params: didn't recognize case. Quitting")
        quit()

   
    
    return mylabel, message, pointsize, pointalpha, mymarker, mylw, chose_colour, specific_colour 






#===================================
def tweak_plot_2D(fig,ax,case):
#===================================
    """
    tweaks the plot, according to case.
    fig:    figure object
    ax:     axis object
    case:   for which case to tweak. Expects a string.


    accepted cases:
        'general':      Nothing particular.


    returns:
        nothing
    """

    #============================
    # Axes, ticks and labelling
    #============================

    if (case=='general'):
        # set tick params (especially digit size)
        ax.tick_params(axis='both', which='major', labelsize=12,top=5)

        #label axes
        ax.set_xlabel(r'x')#, labelpad=20, family='serif',size=16)
        ax.set_ylabel(r'y')#, labelpad=20, family='serif',size=16)


    #=========
    #LEGEND
    #=========

    legloc = 0
    lgnd1=ax.legend(loc=legloc, scatterpoints=1, fancybox=True, framealpha=1)


    #Set dots larger
    for l in range(len(lgnd1.legendHandles)):
        lgnd1.legendHandles[l]._sizes = [30]
        lgnd1.legendHandles[l].set_alpha(1)



    #==================
    # other cosmetics
    #==================
    
    plt.tight_layout()

    return




#=============================================
def tweak_3plot_2D(fig, ax1, ax2, ax3, case):
#=============================================
    """
    tweaks the plot, according to case.
    fig:    figure object
    ax123:     axis object
    case:   for which case to tweak. Expects a string.


    accepted cases:
        'general':      Nothing particular.


    returns:
        nothing
    """

    import numpy as np

    #============================
    # Axes, ticks and labelling
    #============================

    if (case=='general'):
        # set tick params (especially digit size)
        ax1.tick_params(axis='both', which='major', labelsize=12,top=5)
        ax2.tick_params(axis='both', which='major', labelsize=12,top=5)
        ax3.tick_params(axis='both', which='major', labelsize=12,top=5)

        #label axes
        ax1.set_xlabel(r'x')#, labelpad=20, family='serif',size=16)
        ax1.set_ylabel(r'y')#, labelpad=20, family='serif',size=16)

        ax2.set_xlabel(r'y')#, labelpad=20, family='serif',size=16)
        ax2.set_ylabel(r'z')#, labelpad=20, family='serif',size=16)

        ax3.set_xlabel(r'x')#, labelpad=20, family='serif',size=16)
        ax3.set_ylabel(r'z')#, labelpad=20, family='serif',size=16)


    #=========
    #LEGEND
    #=========

    legloc = 0
    lgnd1=ax1.legend(loc=legloc, scatterpoints=1, fancybox=True, framealpha=1)
    lgnd2=ax2.legend(loc=legloc, scatterpoints=1, fancybox=True, framealpha=1)
    lgnd3=ax3.legend(loc=legloc, scatterpoints=1, fancybox=True, framealpha=1)


    #Set dots larger
    for l in range(len(lgnd1.legendHandles)):
        lgnd1.legendHandles[l]._sizes = [30]
        lgnd1.legendHandles[l].set_alpha(1)

    for l in range(len(lgnd2.legendHandles)):
        lgnd2.legendHandles[l]._sizes = [30]
        lgnd2.legendHandles[l].set_alpha(1)

    for l in range(len(lgnd3.legendHandles)):
        lgnd3.legendHandles[l]._sizes = [30]
        lgnd3.legendHandles[l].set_alpha(1)




    #===================================
    # Get all plots to have same size
    #===================================

    xcoords = np.array([])
    ycoords = np.array([])
    zcoords = np.array([])

    nparts = 1
    i = 0

    while True:

        try:
            d1=ax1.collections[i]
            d1.set_offset_position('data')
            xy = d1.get_offsets()

            d2=ax2.collections[i]
            d2.set_offset_position('data')
            yz = d2.get_offsets()
        except IndexError:
            break

        nparts = xy.shape[0]
        if nparts > 0:
            xcoords = np.concatenate((xcoords, xy[:,0]))
            ycoords = np.concatenate((ycoords, xy[:,1]))
            zcoords = np.concatenate((zcoords, yz[:,1]))

        i += 1


    xmin = np.min(xcoords)
    ymin = np.min(ycoords)
    zmin = np.min(zcoords)
    xmax = np.max(xcoords)
    ymax = np.max(ycoords)
    zmax = np.max(zcoords)

    xcenter = 0.5*(xmin + xmax)
    ycenter = 0.5*(ymin + ymax)
    zcenter = 0.5*(zmin + zmax)

    dx = (xmax - xmin)
    dy = (ymax - ymin)
    dz = (zmax - zmin)

    plotsize = max(dx, dy, dz)

    ax1.set_xlim([xcenter-0.5*plotsize, xcenter+0.5*plotsize])
    ax1.set_ylim([ycenter-0.5*plotsize, ycenter+0.5*plotsize])

    ax2.set_xlim([ycenter-0.5*plotsize, ycenter+0.5*plotsize])
    ax2.set_ylim([zcenter-0.5*plotsize, zcenter+0.5*plotsize])

    ax3.set_xlim([xcenter-0.5*plotsize, xcenter+0.5*plotsize])
    ax3.set_ylim([zcenter-0.5*plotsize, zcenter+0.5*plotsize])






    #==================
    # other cosmetics
    #==================
    
    plt.tight_layout()

    return









###########################################################
###########################################################
###########################################################
#################  OLD STUFF  #############################
###########################################################
###########################################################
###########################################################
# 
# 
# #########################
# ##--- PLOTTING 2D      ##
# #########################
# 
# def plot_children2D(children, clumpid, x_part, y_part,z_part,pointsize, pointalpha, toomanychildren, ax1, ax2,ax3):
# 
#     # plot children
#     if (len(children)>0):
#         if (toomanychildren):
#             ax1.scatter(x_part,y_part,s=pointsize,c='blue', label='ptcls of child clumps: '+str(len(x_part)), lw=0, marker=',',alpha=pointalpha)
#             ax2.scatter(y_part,z_part,s=pointsize,c='blue', label='ptcls of child clumps: '+str(len(x_part)), lw=0, marker=',',alpha=pointalpha)
#             ax3.scatter(x_part,z_part,s=pointsize,c='blue', label='ptcls of child clumps: '+str(len(x_part)), lw=0, marker=',',alpha=pointalpha)
# 
# 
#         else:
#                 ax1.scatter(x,y,s=pointsize,c=fullcolorlist[i+1], label='ptcls of child clump '+str(children[i])+':'+str(len(x)), lw=0, marker=',',alpha=pointalpha)
#                 ax2.scatter(y,z,s=pointsize,c=fullcolorlist[i+1], label='ptcls of child clump '+str(children[i])+':'+str(len(x)), lw=0, marker=',',alpha=pointalpha)
#                 ax3.scatter(x,z,s=pointsize,c=fullcolorlist[i+1], label='ptcls of child clump '+str(children[i])+':'+str(len(x)), lw=0, marker=',',alpha=pointalpha)
#     return
# 
# 
# 
# 
# 
# 
# #----------------------------------------------------------
# 
# 
# def get_plot_params(x,which,counter,clumpID):
#     mylw=0
#     mymarker=','
# 
# 
#     if (which=='add-to-halo-cos'):
#         mylabel=None
#         message='halo-namegiver clump'+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=3
#         mylw=0.0
#         pointalpha=0.6
#         mymarker='o'
# 
#     elif (which=='halo'):
#         mylabel='halo-namegiver particles'#+str(clumpID)+': '+str(len(x))
#         message='halo-namegiver clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=.8
#         mylw=0.0
#         pointalpha=0.2
# 
#     if (which=='halo-two'):
#         mylabel='halo-namegiver particles'#+str(clumpID)+': '+str(len(x))
#         message='halo-namegiver clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=1
#         mylw=0.0
#         pointalpha=0.2
#     
#     if (which=='halo-two-phew'):
#         mylabel='halo-namegiver particles'#+str(clumpID)+': '+str(len(x))
#         message='halo-namegiver clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=1
#         mylw=0.0
#         pointalpha=0.2
# 
#     elif (which=='halo-sub'):
#         mylabel='halo-namegiver particles'#+str(clumpID)+': '+str(len(x))
#         message='halo-namegiver clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=.2
#         mylw=0.0
#         pointalpha=0.2
# 
#     elif (which=='halo-cos'):
#         mylabel='halo-namegiver particles'#+str(clumpID)+': '+str(len(x))
#         message='halo-namegiver clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=3
#         mylw=0.0
#         pointalpha=0.6
#         mymarker='o'
# 
#     elif (which=='halo-cos-full'):
#         mylabel='halo particles'#+str(clumpID)+': '+str(len(x))
#         message='halo-namegiver clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=3
#         mylw=0.0
#         pointalpha=0.6
#         mymarker='o'
# 
#     elif (which=='sub-cos'):
#         mylabel='subhalo '+str(counter)+' particles'#+str(clumpid)+':'+str(len(x))
#         message='child clump level'+str(counter )+' : plotting '+str(len(x))+' particles'
#         pointsize=8
#         pointalpha=.9
#         mymarker='o'
#     
# 
#     elif (which=='levels'):
#         mylabel='level '+str(counter-1)+' particles'#+str(clumpID)+':'+str(len(x))
#         message='child clump level'+str(counter -1)+' : plotting '+str(len(x))+' particles'
#         pointsize=1
#         pointalpha=.6
# 
# 
#     elif (which=='sub-two'):
#         mylabel='subhalo particles'#+str(clumpID)+':'+str(len(x))
#         message='child clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=2
#         pointalpha=.8
# 
#     elif (which=='sub-sub'):
#         mylabel='subhalo particles'#+str(clumpID)+':'+str(len(x))
#         message='child clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=0.3
#         pointalpha=.4
# 
#     elif (which=='sub-subsub'):
#         mylabel='subsubhalo '+str(counter )+' particles'#+str(clumpID)+':'+str(len(x))
#         message='child clump '+str(clumpID)+' : plotting '+str(len(x))+' particles'
#         pointsize=.8
#         pointalpha=.8
# 
#    
#     
#     return mylabel, message, pointsize,pointalpha, mymarker, mylw
#    
#     
# 
# 
# 
# 
# def plot_2D(x,y,ax1,which,counter):
#   
# 
#     mylabel, message, pointsize,pointalpha,mymarker,mylw=get_plot_params(x,which,counter,-1)
# 
#     print mylabel, message, pointsize, pointalpha, mylabel, mylw
#     print message
#     ax1.scatter(x,y,s=pointsize,c=fullcolorlist[counter], label=mylabel, lw=mylw, marker=mymarker,alpha=pointalpha) 
#     # ax1.scatter(x,y) 
# 
#     return
# 
# 
# 
# 
# 
# 
# #----------------------------------------------------------
# 
# 
# 
# 
# 
# 
# 
# def make_fancy_axes(ax1,ax2,ax3,ax4,clumpx, clumpy, clumpz, x_part, y_part, z_part,twodim):
#     #needs clumpx, clumpy,clumpz from get_clump_data()
#     #gets corrections for various axes
#     print "Calculating plotting domain."
#     maxx=0.0
#     maxy=0.0
#     maxz=0.0
#     for i in range(len(x_part)):
#         distx=abs(x_part[i]-clumpx)
#         if (distx>maxx):
#             maxx=distx
#         disty=abs(y_part[i]-clumpy)
#         if (disty>maxy):
#             maxy=disty
#         distz=abs(z_part[i]-clumpz)
#         if (distz>maxz):
#             maxz=distz
# 
# 
# 
# 
# 
#     if (twodim):
#         if(maxx>=maxy):
#             xyc=maxx
#         else:
#             xyc=maxy
#         
#         if(maxy>=maxz):
#             yzc=maxy
#         else:
#             yzc=maxz
#         
#         if(maxx>=maxz):
#             xzc=maxx
#         else:
#             xzc=maxy
# 
#         if (xyc>0):
#             ax1.set_xlim(clumpx-1.5*xyc, clumpx+1.5*xyc)
#             ax1.set_ylim(clumpy-1.5*xyc, clumpy+1.5*xyc)   
#         if (yzc>0):
#             ax2.set_xlim(clumpy-1.5*yzc, clumpy+1.5*yzc)
#             ax2.set_ylim(clumpz-1.5*yzc, clumpz+1.5*yzc)   
#         if (xzc>0):
#             ax3.set_xlim(clumpx-1.5*xzc, clumpx+1.5*xzc)
#             ax3.set_ylim(clumpz-1.5*xzc, clumpz+1.5*xzc)   
# 
#     else: # 3-dim
#         corr=max((maxx,maxy,maxz))
#         if (corr>0):
#             ax1.set_xlim(clumpx-1.5*corr, clumpx+1.5*corr)
#             ax1.set_ylim(clumpy-1.5*corr, clumpy+1.5*corr)
#             ax1.set_zlim(clumpz-1.5*corr, clumpz+1.5*corr)
#             ax2.set_xlim(clumpx-1.5*corr, clumpx+1.5*corr)
#             ax2.set_ylim(clumpy-1.5*corr, clumpy+1.5*corr)
#             ax2.set_zlim(clumpz-1.5*corr, clumpz+1.5*corr)
#             ax3.set_xlim(clumpx-1.5*corr, clumpx+1.5*corr)
#             ax3.set_ylim(clumpy-1.5*corr, clumpy+1.5*corr)
#             ax3.set_zlim(clumpz-1.5*corr, clumpz+1.5*corr)
#             ax4.set_xlim(clumpx-1.5*corr, clumpx+1.5*corr)
#             ax4.set_ylim(clumpy-1.5*corr, clumpy+1.5*corr)
#             ax4.set_zlim(clumpz-1.5*corr, clumpz+1.5*corr)
#     return 
# 
# 
# 
# 
# 
# 
# #----------------------------------------------------------
# 
# 
# 
# 
# 
# 
# def tweak_plot_2D(fig,plt,ax1,case):
#     import matplotlib.pyplot as plt
#     from matplotlib.font_manager import FontProperties # for legend
#      # print ax1.azim, ax1.elev, ax1.dist
# 
#     fontP=FontProperties()
#     fontP.set_size('x-large') 
#     fontP.set_family('serif') # families = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace']
# 
#     #shift axes
#     if (case=='cos'):
#         # set tick params (especially digit size)
#         ax1.tick_params(axis='both',which='major',labelsize=12,top=5)
# 
#         #label axes
#         ax1.set_xlabel(r'x', labelpad=15, family='serif',size=16)
#         ax1.set_ylabel(r'y', labelpad=15, family='serif',size=16)
#         # ax1.set_zlabel(r'kpc', labelpad=15, family='serif',size=16)
#         ax1.set_xlim(0,1)
#         ax1.set_ylim(0,1)
#         plt.subplots_adjust(left=0., right=1, top=1, bottom=0.0,wspace=0.0,hspace=0.0)
# 
#     if (case=='mf'):
#         ax1.tick_params(axis='both',which='major',labelsize=12,top=5)
#         ax1.set_xlabel(r'subclump particles (binned)', labelpad=15, family='serif',size=16)
#         ax1.set_ylabel(r'number of subclumps in bin', labelpad=15, family='serif',size=16)
#     
#         plt.subplots_adjust(left=0.5, right=.95, top=.95, bottom=0.5,wspace=0.0,hspace=0.0)
# 
#         fontP.set_size('medium') 
#     
#     
# 
# 
# 
#     #LEGEND
# 
# 
#     lgnd1=ax1.legend(loc=0, scatterpoints=1,prop=fontP, framealpha=1)
#     for l in range(len(lgnd1.legendHandles)):
#         lgnd1.legendHandles[l]._sizes = [20]
#         lgnd1.legendHandles[l].set_alpha(1)
# 
# 
# 
# 
# 
# 
# 

