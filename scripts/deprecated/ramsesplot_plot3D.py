
# A module for the ramsesplot scripts
# Containts all functions about 3D plots

# Contains:
#   setup_fig
#   save_fig
#   plot_3D
#   get_plot_params
#   tweak_plot_3D
# 

import matplotlib.pyplot as plt
from matplotlib.pyplot import subplots_adjust
from matplotlib.font_manager import FontProperties # for legend
from mpl_toolkits.mplot3d import Axes3D



##########################
##--- PLOTTING          ##
##########################




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
    ax = fig.add_subplot(111, projection='3d')


    return fig, ax 










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
def plot_3D(x, y, z, ax, case, clumpID):
#============================================
    """
    plot given x,y,z. First gets parameters for case.
    x, y, z:    x, y, z coordinates to plot
    ax :        axis object
    case:       for which case to plot/choose parameters
    clumpID:    clump ID which is to be plotted.


    accepted cases:
        'halo'      halo namegiver particles
        'child'     child clump particles

    returns:
        nothing

    """

    mylabel, message, pointsize, pointalpha, mymarker, mylw = get_plot_params(case, x.shape[0], clumpID)

    print(message)
    ax.scatter(x, y, z, 
            s=pointsize, 
            label=mylabel, 
            lw=mylw, 
            marker=mymarker, 
            depthshade=True, 
            alpha=pointalpha)
    return








#============================================
def get_plot_params(case, length, clumpID):
#============================================
    """
    set up plot parameters for possible cases.
    case:      which case to plot.
               passed on by the plot_3D function.
    length:    lenght of array to plot
    clumpID:   clump ID

    returns:
        mylabel, message, pointsize, pointalpha, mymarker, mylw

    mylabel:    legend label for plot
    message:    message to print before plotting
    pointsize:  size of plot points
    pointalpha: alpha of points to plot
    mymarker:   plot markers
    mylw:       line width around plot points

    """



    #default values

    mylabel=None
    message='No message'
    pointsize=1
    mylw=0.0
    pointalpha=0.8
    mymarker=','


    if (case=='child'):
        mylabel='child ' + str(clumpID)
        message=mylabel +' : plotting ' + str(length) + ' particles'
        pointsize=4
        mylw=0.0
        pointalpha=0.1
        mymarker='o'


    elif (case=='halo'):
        mylabel= 'namegiver '+str(clumpID)
        message='halo ' + mylabel + ' : plotting ' +str(length) + ' particles'
        pointsize=1
        mylw=0.0
        pointalpha=0.1
        mymarker=','

    else:
        print("get_plot_params: didn't recognize case. Quitting")
        quit()

   
    
    return mylabel, message, pointsize, pointalpha, mymarker, mylw









#===================================
def tweak_plot_3D(fig,ax,case):
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
        ax.set_zlabel(r'z')#, labelpad=20, family='serif',size=16)


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



###########################################################
###########################################################
###########################################################
#################  OLD STUFF  #############################
###########################################################
###########################################################
###########################################################



# #########################
# ##--- PLOTTING 3D      ##
# #########################
# 
# 
# def plot_3D(x,y,z,ax1,which,counter,clumpID):
# 
#     from mpl_toolkits.mplot3d import Axes3D
# 
#     mylabel, message, pointsize,pointalpha, mymarker, mylw=get_plot_params(x,which,counter,-1)
# 
#     print message
#     ax1.scatter(x,y,z,s=pointsize,c=fullcolorlist[counter], label=mylabel, lw=mylw, marker=mymarker,depthshade=True,alpha=pointalpha)
#     return
# 
# 
# #----------------------------------------------------------
# 
# 
# 
# 
# 
# 
# def tweak_plot_3D(fig,plt, ax1,case):
#     import matplotlib.pyplot as plt
#     from matplotlib.pyplot import subplots_adjust
#     from matplotlib.font_manager import FontProperties # for legend
#      # print ax1.azim, ax1.elev, ax1.dist
# 
# 
#     # default Legend parameters
#     legloc = 0
#     legsize = 'x-large'
#     bbox=False
# 
#     #shift axes
#     if (case=='cos'): 
#         ax1.view_init(elev=40,azim=40)
#         ax1.tick_params(axis='both',which='major',labelsize=18,pad=12)
# 
#         ax1.set_xlabel(r'x', labelpad=20, family='serif',size=24)
#         ax1.set_ylabel(r'y', labelpad=20, family='serif',size=24)
#         ax1.set_zlabel(r'z', labelpad=25, family='serif',size=24)
# 
# 
#         # ax1.set_xticklabels(["",-0.015,"","","","","",0.015,""],va='baseline',ha='right')
#         # ax1.set_yticklabels(["","0.310","","0.330",""],va='baseline',ha='left')
#         # ax1.set_zticklabels(["","0.885","","","","","","0.915",""],va='baseline',ha='left')
# 
#         ax1.set_xlim3d(-0.02,0.02)
#         ax1.set_ylim3d(0.30,0.35)
#         ax1.set_zlim3d(0.880,0.920)
#         
#         
#         ax1.w_xaxis.set_pane_color( (0,0.1,0.7,.15) )
#         ax1.w_yaxis.set_pane_color( (0,0.1,0.7,.15) )
#         ax1.w_zaxis.set_pane_color( (0,0.1,0.7,.15) )
# 
#         ax1.grid()
#         ax1.w_xaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_yaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_zaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
# 
#         legloc = 'upper left'
#         legsize = 22
#         bbox=True
#         bbox_to_anchor=(.51, 1.01)
# 
#         subplots_adjust(left=-0.01, right=1.01, top=1.00, bottom=-0.01,wspace=0.00,hspace=0.00)
# 
# 
# 
# 
#     elif (case=='two'):
#         # ax1.view_init(elev=5,azim=-110)
#         ax1.view_init(elev=30,azim=290)
# 
#         # set tick params (especially digit size)
#         ax1.tick_params(axis='both',which='major',labelsize=18,top=5)
# 
#         #label axes
#         ax1.set_xlabel(r'x $[kpc]$', labelpad=5, family='serif',size=24)
#         ax1.set_ylabel(r'y $[kpc]$', labelpad=5, family='serif',size=24)
#         ax1.set_zlabel(r'z $[kpc]$', labelpad=5, family='serif',size=24)
#         
#         ax1.set_xticklabels([0,"","","","",500],va='baseline',ha='right')
#         ax1.set_yticklabels([0,"","","","",500],va='baseline',ha='left')
#         ax1.set_zticklabels([0,"","","","",500],va='baseline',ha='left')
#         ax1.set_xlim3d(0,500)
#         ax1.set_ylim3d(0,500)
#         ax1.set_zlim3d(0,500)
# 
#         ax1.grid()
#         ax1.w_xaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_yaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_zaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_xaxis.set_pane_color(( 0,0.1,0.7,.15))
#         ax1.w_yaxis.set_pane_color(( 0,0.1,0.7,.15))
#         ax1.w_zaxis.set_pane_color(( 0,0.1,0.7,.15))
# 
# 
#         legloc = 'upper left'
#         # legsize = 'xx-large'
#         legsize = 22
# 
#         bbox=True
#         subplots_adjust(left=-0.10, right=1.03, top=1.03, bottom=-0.03,wspace=0.00,hspace=0.00)
#         # bbox_to_anchor=(0.55, .97)
#         bbox_to_anchor=(0.55, 0.98)
# 
# 
# 
#     elif (case=='sub'):
#         # ax1.view_init(azim=5, elev=15.)
#         # ax1.view_init(elev=60,azim=340)
#         ax1.view_init(azim=310, elev=160.)
# 
#         ax1.set_xlim3d(0,800)
#         ax1.set_ylim3d(0,800)
#         ax1.set_zlim3d(0,800)
# 
#         ax1.grid(True)
#         ax1.w_xaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_yaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_zaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_xaxis.set_pane_color(( 0,0.1,0.7,.15))
#         ax1.w_yaxis.set_pane_color(( 0,0.1,0.7,.15))
#         ax1.w_zaxis.set_pane_color(( 0,0.1,0.7,.15))
# 
# 
# 
# 
#         # set tick params (especially digit size)
#         ax1.tick_params(axis='both',which='major',labelsize=18,pad=0)
# 
#         #label axes
#         ax1.set_xlabel(r'x $[kpc]$', labelpad=2, family='serif',size=24)
#         ax1.set_ylabel(r'y $[kpc]$', labelpad=2, family='serif',size=24)
#         # ax1.zaxis.set_rotate_label(False)
#         ax1.set_zlabel(r'z $[kpc]$', labelpad=5, family='serif',size=24)#,rotation=90)
# 
#         ax1.set_xticklabels([0,"","","","","","","",800],va='bottom',ha='left')
#         ax1.set_yticklabels([0,"","","","","","","",800],va='baseline',ha='right')
#         ax1.set_zticklabels([0,"","","","","","","",800],va='baseline',ha='right')
# 
# 
# 
#         legloc = 'upper left'
#         legsize = 22
# 
#         subplots_adjust(left=-0.03, right=1.03, top=1.02, bottom=0.00,wspace=0.00,hspace=0.00)
# 
#         bbox=True
#         bbox_to_anchor=(0.52, .23)
# 
# 
# 
#     elif (case=='sub-sub'):
# 
# 
#         # ax1.view_init(azim=5, elev=15.)
#         # ax1.view_init(elev=60,azim=340)
#         ax1.view_init(azim=310, elev=160.)
# 
#         ax1.set_xlim3d(300,800)
#         ax1.set_ylim3d(300,800)
#         ax1.set_zlim3d(300,800)
# 
#         ax1.grid(True)
#         ax1.w_xaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_yaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_zaxis._axinfo.update({'grid' : {'color': (1, 1, 1, 1)}})
#         ax1.w_xaxis.set_pane_color(( 0,0.1,0.7,.15))
#         ax1.w_yaxis.set_pane_color(( 0,0.1,0.7,.15))
#         ax1.w_zaxis.set_pane_color(( 0,0.1,0.7,.15))
# 
# 
# 
# 
#         # set tick params (especially digit size)
#         ax1.tick_params(axis='both',which='major',labelsize=18,pad=0)
# 
#         #label axes
#         ax1.set_xlabel(r'x $[kpc]$', labelpad=2, family='serif',size=24)
#         ax1.set_ylabel(r'y $[kpc]$', labelpad=2, family='serif',size=24)
#         # ax1.zaxis.set_rotate_label(False)
#         ax1.set_zlabel(r'z $[kpc]$', labelpad=5, family='serif',size=24)#,rotation=90)
# 
#         ax1.set_xticklabels([300,"","","","",800],va='bottom',ha='left')
#         ax1.set_yticklabels([300,"","","","",800],va='baseline',ha='right')
#         ax1.set_zticklabels([300,"","","","",800],va='baseline',ha='right')
# 
# 
# 
#         legloc = 'upper left'
#         legsize = 22
# 
#         subplots_adjust(left=-0.03, right=1.03, top=1.02, bottom=0.00,wspace=0.00,hspace=0.00)
# 
#         bbox=True
#         bbox_to_anchor=(0.52, .23)
# 
# 
# 
#   
# 
# 
# 
#     #LEGEND
# 
#     fontP=FontProperties()
#     fontP.set_size(legsize)
#     fontP.set_family('serif') # families = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace']
#     if bbox:
#         lgnd1=ax1.legend(scatterpoints=1,prop=fontP, fancybox=True, framealpha=1, bbox_to_anchor=bbox_to_anchor)
#     else:
#         lgnd1=ax1.legend(loc=legloc, scatterpoints=1,prop=fontP, fancybox=True, framealpha=1)
# 
#     for l in range(len(lgnd1.legendHandles)):
#         lgnd1.legendHandles[l]._sizes = [30]
#         lgnd1.legendHandles[l].set_alpha(1)
# 
#        
# 
# 
# def save_fig(this_name,fig,workdir):
#     import matplotlib.pyplot as plt
#     fig_path = workdir+'/'+this_name+'.png'
#     print "saving figure as "+fig_path
#     plt.savefig(fig_path, format='png', facecolor=fig.get_facecolor(), transparent=False, dpi=100)#,bbox_inches='tight' )
#     print "done", this_name+'.png'
# 
# 
# 
# 
# 
# #----------------------------------------------------------
# #----------------------------------------------------------
# #----------------------------------------------------------
# 









