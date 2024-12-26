# A module for the ramsesplot scripts
# Containts all functions about getting data

# Contains:
#   get_clump_data
#   get_children_data
#   get_level_data
#   get_halo_particles
#   get_all_halo_particles
#   get_all_level_particles

# functions not to be used externally:
#   read_particles

import numpy as np
from sys import argv #command line arguments
import warnings


#########################
##--- ACQUIRING DATA   ##
#########################



#============================================================
def get_all_halo_particles(noutput,particles):
#============================================================
    """
    reads through particle data and extracts all particles that are in a
    halo/clump.

    noutput:    number of files to go through
    particles:  number of particles of sim in total


    returns:
        hx, hy, hz

    hx, hy, hz: numpy array of halo particles.

    """

    # read in particles
    data = read_particles(3, particles, noutput)

    # filter out all non-clump particles
    # store child particles in array for each child.
    # chose all particles for first guess of length.
    halox=np.zeros(particles)
    haloy=np.zeros(particles)
    haloz=np.zeros(particles)


    # set up found particle counter
    hc=0

    for i in range(data.shape[0]):
        if (data[i,3] != 0):
            #find halo particles
            # append halo particles 
            halox[hc]=data[i,0]
            haloy[hc]=data[i,1]
            haloz[hc]=data[i,2]
            hc+=1
             

    hx=np.zeros(hc)
    hy=np.zeros(hc)
    hz=np.zeros(hc)
    for i in range(hc):
        hx[i]=halox[i]
        hy[i]=haloy[i]
        hz[i]=haloz[i]


    return hx, hy, hz




#======================================================================
def get_all_level_particles(noutput, particles, clump_levels, halos):
#======================================================================
    """
    Sorts all particles by level.

    noutput:        number of files to go through
    particles:      number of particles of sim in total
    clump_levels:   list of all clumps; clump_levels[2] = all level 2 clumps
    halos:          list of all halos


    returns:
        px, py, pz, hx, hy, hz

    px, py, pz: list of numpy arrays of satellite halo particles, aranged by level
    hx, hy, hz: list of numpy arrays numpy array of halo particles, aranged by level

    """

    # read in particles
    data = read_particles(3+noutput, particles, noutput)



    # filter out all non-clump particles
    # store child particles in array for each child.
    # chose all particles for first guess of length.
    
    maxlevel = len(clump_levels)
    halox=np.zeros(particles)
    haloy=np.zeros(particles)
    haloz=np.zeros(particles)
    

    parx = []
    pary = []
    parz = []
    for l in range(maxlevel):
        parx.append(np.zeros(particles))
        pary.append(np.zeros(particles))
        parz.append(np.zeros(particles))


    # set up found particle counter
    hc=0
    pc=[0]*maxlevel




    print("Sorting particles by levels.")

    for i in range(data.shape[0]):
        clumpID = data[i,3]
        if (clumpID != 0): #particle is in clump
            if clumpID in halos:
                halox[hc]=data[i,0]
                haloy[hc]=data[i,1]
                haloz[hc]=data[i,2]
                hc+=1
            else:
                for level in range(maxlevel):
                    if clumpID in clump_levels[level]:
                        parx[level][pc[level]]=data[i,0]
                        pary[level][pc[level]]=data[i,1]
                        parz[level][pc[level]]=data[i,2]
                        pc[level]+=1
                        break
                    

             

    hx=np.zeros(hc)
    hy=np.zeros(hc)
    hz=np.zeros(hc)
    for i in range(hc):
        hx[i]=halox[i]
        hy[i]=haloy[i]
        hz[i]=haloz[i]
        
    
    px = []
    py = []
    pz = []
    
    for l in range(maxlevel):
        px.append(np.zeros(pc[l]))
        py.append(np.zeros(pc[l]))
        pz.append(np.zeros(pc[l]))
        
        for p in range(pc[l]):
            px[l][p] = parx[l][p]
            py[l][p] = pary[l][p]
            pz[l][p] = parz[l][p]                   
        
    return px, py, pz, hx, hy, hz




#===================================
def get_children_data(halo,noutput):
#===================================
    """
    Finds children and their levels for given halo. 

    halo:       halo ID for which to look for
    noutput:    number of files through which to look through

    returns:
        children, child_levels
        children:       sorted list of all children of halo
        child_levels:   list of levels of nchildren, sorted appropriately
    """

    # Get list of children for halo
    # noutput: number of outputfiles


    # read in clump data
    data_int, data = get_clump_data(4,noutput)
    
    
    #identify all children and their levels recursively for given halo
    children=[]
    child_levels=[]
    i=0
    while(i<data_int.shape[0]):
        try:
            if (data_int[i,2] in children or data_int[i,2]==halo): #if parent is halo or one of the children
                if (data_int[i,0] != halo): #if clump itself isn't halo
                    if (data_int[i,0] not in children):
                        children.append(data_int[i,0])
                        child_levels.append(data_int[i,1])
                        i=-1                #start again from the beginning to find hierarchical buildup
        
        except IndexError: # in case there is only one clump in simulation
            children = []
            child_levels = []
            break

        i+=1

    # sort children list and child_levels list appropriately
    # first create a list of tuples
    together = zip(children, child_levels)
    # then extract correct value for the sorted lists.
    # tuples are sorted by the first value.
    children_sorted = [x[0] for x in sorted(together)]
    child_levels_sorted = [x[1] for x in sorted(together)]

    # return clumpx, clumpy, clumpz, children
    return children_sorted, child_levels_sorted




#======================================
def get_clump_data(start_ind, noutput):
#======================================
    """
    reads in clumpfinder data from files.
    start_ind:      index from which to start out from when reading filenames 
                    given as cmd line args
    noutput:        number of files to go through
    
    returns:
        data:     numpy array of the data
            data[0]: x coordinate of clump
            data[1]: y coordinate of clump
            data[2]: z coordinate of clump
        data_int: numpy array of integer data
            data_int[0]: clump ID
            data_int[1]: clump level
            data_int[2]: ID of parent clump
    """


    print("Reading in clump data.")


    # Loop over files
    for i in range(0,noutput):
        inputfile=str(argv[i+start_ind]) 
        
        # get clump center
        # ignore empty files
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            temp_data=np.loadtxt(inputfile, dtype='float', skiprows=1, usecols=[4,5,6])
        if(temp_data.shape[0]>0):
            if 'data' in locals():
                data = np.vstack((data, temp_data))
            else:
                data = temp_data
        
        #get clump ids, parents and levels
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            temp_data_int=np.loadtxt(inputfile, dtype='int', skiprows=1, usecols=[0,1,2])
        if(temp_data_int.shape[0]>0):
            if 'data_int' in locals():
                data_int= np.vstack((data_int, temp_data_int))
            else:
                data_int = temp_data_int


    return data_int, data




#=======================================================
def get_halo_particles(children,halo,noutput,particles):
#=======================================================
    """
    reads through particle data and extracts the ones for the halo.
    children:   list of children
    halo:       halo ID
    noutput:    number of files to go through
    particles:  number of particles of sim in total


    returns:
        fx, fy, fz, hx, hy, hz

    fx, fy, fz: lists of numpy arrays of particles. 
                fx[i] is the particle list for child children[i]
    hx, hy, hz: numpy array of halo particles.

    """


    # read in particles
    data = read_particles(4+noutput, particles, noutput)

    nchildren = len(children)
    
    
    # filter out all non-clump particles
    # store child particles in array for each child.
    # chose all particles for first guess of length.
    filteredx = []
    filteredy = []
    filteredz = []

    for c in range(nchildren):
        filteredx.append(np.zeros(particles))
        filteredy.append(np.zeros(particles))
        filteredz.append(np.zeros(particles))

    halox=np.zeros(particles)
    haloy=np.zeros(particles)
    haloz=np.zeros(particles)


    # set up found particle counters
    fc=nchildren*[0]
    hc=0

    for i in range(particles):
        if (len(children)>0):
        #find children particles
            if abs(data[i,3]) in children:
                childind = children.index(abs(data[i,3]))
                # if particle belongs to child, append to correct array
                filteredx[childind][fc[childind]]=data[i,0]  
                filteredy[childind][fc[childind]]=data[i,1]
                filteredz[childind][fc[childind]]=data[i,2] 
                # raise found particle index
                fc[childind]+=1


        if (abs(data[i,3])==halo):
            #find halo particles
            # append halo particles 
            halox[hc]=data[i,0]
            haloy[hc]=data[i,1]
            haloz[hc]=data[i,2]
            hc+=1
            

   

    #create empty lists for final particle lists
    fx = []
    fy = []
    fz = []

    #add array for each child with correct number of particles
    for c in range(nchildren):
        fx.append(np.zeros(fc[c]))
        fy.append(np.zeros(fc[c]))
        fz.append(np.zeros(fc[c]))

    #copy the particles into final arrays
    for c in range(nchildren):
        for i in range(fc[c]):
            fx[c][i]=filteredx[c][i]
            fy[c][i]=filteredy[c][i]
            fz[c][i]=filteredz[c][i]

    hx=np.zeros(hc)
    hy=np.zeros(hc)
    hz=np.zeros(hc)
    for i in range(hc):
        hx[i]=halox[i]
        hy[i]=haloy[i]
        hz[i]=haloz[i]


    return fx, fy, fz, hx, hy, hz




#===================================
def get_level_data(noutput):
#===================================
    """
    Finds all clumps, gives them sorted by level and halos.

    noutput:    number of files through which to look through

    returns:
        clump_levels, halos
        clump_levels:      list of all satellite haloes, sorted by level:
                           clump_levels[2] gives list of all level 2 satellites
        halos:             list of all main halos
    """

    # Get list of children for halo
    # noutput: number of outputfiles

 
    # read in clump data
    data_int, data = get_clump_data(3, noutput)
    
    # sort clumps by level or to halo list if halo
    maxlevel = max(data_int[:,1])

    clump_level = []
    for i in range(maxlevel + 1):
        clump_level.append([])
    halos = []
    
    for i in range(data_int.shape[0]):
        #check if halo
        if data_int[i,0]==data_int[i,2]:
            halos.append(data_int[i,0])
        else:
            level = data_int[i,1]
            clump_level[level].append(data_int[i,0])
    
    
    return clump_level, halos




#============================================================
def read_particles(start_ind, particles, noutput):
#============================================================
    """
    extracts all particle data from files. 
    start_ind: index from which to start out from when reading filenames given as cmd line args
    particles:  number of particles of sim in total
    noutput:        number of files to go through
    
    returns:
        data:     numpy array of the data
            data[0]: x coordinate of particles
            data[1]: y coordinate of particles
            data[2]: z coordinate of particles
            data[3]: clump ID the particle is assigned to

    """


    print("Reading in particle data.")
    
    data=np.zeros((particles,4))
    particle_index=0

    for i in range(0,noutput):
        inputfile=str(argv[i+start_ind]) 
        temp_data=np.loadtxt(inputfile, dtype='float', skiprows=1, usecols=[0,1,2,6])
    
        if (temp_data.shape[0]>0):
            for k in range(4):
                for j in range(temp_data.shape[0]):
                    data[j+particle_index,k]=temp_data[j,k]
        
            particle_index += temp_data.shape[0]-1


    return data






###########################################################
###########################################################
###########################################################
#################  OLD STUFF  #############################
###########################################################
###########################################################
###########################################################







##########################
##--- GENERAL VARIABLES ##
##########################


# fullcolorlist=['red', 'green', 'blue', 'black', 'magenta', 'lime','cyan','mediumpurple', 'gold','lightpink','lime','saddlebrown','darkolivegreen','cornflowerblue','dimgrey','navajowhite','black','darkslategray','silver','mediumseagreen','orange','midnightblue','silver']












#########################
##--- ACQUIRING DATA   ##
#########################


# def get_level_list(noutput):
#     for i in range(0,noutput):
#         inputfile=str(argv[i+4]) 
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             temp_data_int=np.loadtxt(inputfile, dtype='int', skiprows=1, usecols=[0,1])
#         if(temp_data_int.shape[0]>0):
#             if 'data_int' in locals():
#                 data_int= np.vstack((data_int, temp_data_int))
#             else:
#                 data_int = temp_data_int
# 
#     
#     maxid_clump=max(data_int[:,0]) #how many clumps
# 
#     for i in range(0,noutput):
#         inputfile=str(argv[i+4+2*noutput]) 
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             temp_data_int=np.loadtxt(inputfile, dtype='int', skiprows=1, usecols=[0])
#         if(len(np.atleast_1d(temp_data_int))>1):# and temp_data_int.shape[0]>0):
#             if 'halo_data' in locals():
#                 halo_data= np.concatenate((halo_data, temp_data_int))
#             else:
#                 halo_data = temp_data_int
# 
#     maxid_halo=max(halo_data) #how many clumps
#     clumpnr=max((maxid_halo, maxid_clump))
# 
#     levels=np.zeros(clumpnr+1)
#     for k in range(len(levels)):
#         levels[k]=-2
#     
# 
#     for j in range(data_int.shape[0]):
#         # if data_int[j][0]==data_int[j][2]:
#         #     #is halo: set -1
#         #     levels[data_int[j][0]]=-1
#         # else:
#             # is child: set clump level
#         levels[data_int[j][0]]=data_int[j][1]
# 
#     for k in range(len(halo_data)):
#         levels[int(halo_data[k])]=-1
# 
#     maxlevel=int(max(levels))
#     # return halos,maxlevel,children,child_levels
#     return levels, maxlevel
# 
# 
# 
# 
# #----------------------------------------------------------
# 
# def get_particle_data_level(levels,noutput,particles,maxlevel):
#     # for cosmo scripts; find particles of halos from halo list,
#     # store them by level
#     # READ CLUMP DATA IN FIRST
#     print "Reading in particle data"
#     data=np.zeros((particles,5))
#     particle_index=0
#     for i in range(0,noutput):
#         inputfile=str(argv[i+4+noutput]) 
#         temp_data=np.loadtxt(inputfile, dtype='float', skiprows=1, usecols=[0,1,2,6,7])
#     
#         if (temp_data.shape[0]>0):
#             for j in range(temp_data.shape[0]):
#                 for k in range(5):
#                     data[j+particle_index,k]=temp_data[j,k]
#         
#             particle_index=particle_index+temp_data.shape[0]
# 
#     c=0
#     for j in range(data.shape[0]):
#         if data[j][3]==0:
#             c+=1
#     print "counting nonbelonging particles ", c
# 
#     print data.shape
#     print particle_index
#     levelind=np.zeros(maxlevel+1)
#     x=[]
#     y=[]
#     z=[]
#     for i in range(maxlevel+1):
#         x.append(np.zeros(particles))
#         y.append(np.zeros(particles))
#         z.append(np.zeros(particles))
#     
#     halox=np.zeros(particles)
#     haloy=np.zeros(particles)
#     haloz=np.zeros(particles)
#     hc=0
# 
# 
# 
#     for i in range(data.shape[0]):
#         if levels[int(data[i,3])]>=0 :
#             lev=int(levels[int(data[i,3])])
#             ind=int(levelind[lev])
#             x[lev][ind]=(data[i,0])
#             y[lev][ind]=(data[i,1])
#             z[lev][ind]=(data[i,2])
#             levelind[lev]+=1
# 
# 
#         elif (levels[int(data[i,3])]==-1):
#             # append halo particles 
#             halox[hc]=(data[i,0])
#             haloy[hc]=(data[i,1])
#             haloz[hc]=(data[i,2])
#             hc+=1
# 
# 
#     print levelind
#     print hc
# 
#     x_final=[]
#     y_final=[]
#     z_final=[]
#     for i in range(maxlevel+1):
#         x_final.append(np.zeros(int(levelind[i])))
#         y_final.append(np.zeros(int(levelind[i])))
#         z_final.append(np.zeros(int(levelind[i])))
#         for j in range(int(levelind[i])):
#             x_final[i][j]=x[i][j]
#             y_final[i][j]=y[i][j]
#             z_final[i][j]=z[i][j]
# 
#     halox_final=np.zeros(hc) 
#     haloy_final=np.zeros(hc) 
#     haloz_final=np.zeros(hc) 
# 
#     
#     for i in range(hc):
#         halox_final[i]=halox[i]
#         haloy_final[i]=haloy[i]
#         haloz_final[i]=haloz[i]
# 
#         
# 
# 
#     return x_final,y_final,z_final,halox_final,haloy_final,haloz_final
# 
# 
# 
# 
# #==========================================
# #==========================================
# #==========================================
# 
# def get_particle_data(children,halo,noutput,particles):
#     # READ CLUMP DATA IN FIRST
#     print "Reading in particle data"
#     data=np.zeros((particles,5))
#     particle_index=0
#     for i in range(0,noutput):
#         inputfile=str(argv[i+4+noutput]) 
#         temp_data=np.loadtxt(inputfile, dtype='float', skiprows=1, usecols=[0,1,2,6,7])
#     
#         if (temp_data.shape[0]>0):
#             for j in range(temp_data.shape[0]):
#                 for k in range(5):
#                     data[j+particle_index,k]=temp_data[j,k]
#         
#             particle_index=particle_index+temp_data.shape[0]-1
# 
#     # filter out all non-clump particles
#     filteredx=np.zeros(particles) #first guess for length.
#     filteredy=np.zeros(particles)
#     filteredz=np.zeros(particles)
#     clumpid=np.zeros(particles)
#     halox=np.zeros(particles)
#     haloy=np.zeros(particles)
#     haloz=np.zeros(particles)
#     unboundx=np.zeros(particles)
#     unboundy=np.zeros(particles)
#     unboundz=np.zeros(particles)
#     unboundclumpid=np.zeros(particles)
# 
#     fc=0
#     hc=0
#     uc=0
# 
#     for i in range(data.shape[0]):
#         if (len(children)>0):
#             # for j in range(0,len(children)):
#                 # if (data[i,3]==children[j]):
#             if (data[i,3] in children):
#                 filteredx[fc]=data[i,0]  
#                 filteredy[fc]=data[i,1]
#                 filteredz[fc]=data[i,2] 
#                 clumpid[fc]=data[i,3]
#                 fc+=1
# 
#                 #get unbound prtcls
#                 if (data[i,4]>0):
#                     unboundx[uc]=data[i,0]
#                     unboundy[uc]=data[i,1]
#                     unboundz[uc]=data[i,2]
#                     unboundclumpid[uc]=data[i,3]
#                     uc+=1
# 
#         if (data[i,3]==halo):
#             # append halo particles 
#             halox[hc]=data[i,0]
#             haloy[hc]=data[i,1]
#             haloz[hc]=data[i,2]
#             hc+=1
#             
#             #get unbound prtcls
#             if (data[i,4]>0):
#                 unboundx[uc]=data[i,0]
#                 unboundy[uc]=data[i,1]
#                 unboundz[uc]=data[i,2]
#                 unboundclumpid[uc]=data[i,3]
#                 uc+=1
# 
# 
#     fx=np.zeros(fc)
#     fy=np.zeros(fc)
#     fz=np.zeros(fc)
#     cid=np.zeros(fc)
#     for i in range(fc):
#         fx[i]=filteredx[i]
#         fy[i]=filteredy[i]
#         fz[i]=filteredz[i]
#         cid[i]=clumpid[i]
# 
#     hx=np.zeros(hc)
#     hy=np.zeros(hc)
#     hz=np.zeros(hc)
#     for i in range(hc):
#         hx[i]=halox[i]
#         hy[i]=haloy[i]
#         hz[i]=haloz[i]
# 
#     ux=np.zeros(uc)
#     uy=np.zeros(uc)
#     uz=np.zeros(uc)
#     ucid=np.zeros(uc)
#     for i in range(uc):
#         ux[i]=unboundx[i]
#         uy[i]=unboundy[i]
#         uz[i]=unboundz[i]
#         ucid[i]=unboundclumpid[i]
# 
#     return fx, fy, fz, cid, hx, hy, hz, ux, uy, uz, ucid
# 
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
# def get_level_particles(child_levels,children,ilevel,x_part,y_part,z_part,clumpid):
#     x=None
#     y=None
#     z=None
#    
#     for i in range(len(children)):
#         if (child_levels[i]==ilevel):
#             xch,ych,zch=get_child_particles(x_part,y_part,z_part,clumpid,children[i])
#             if (xch.shape[0]>0):
#                 if x is None:
#                     x=xch
#                     y=ych
#                     z=zch
#                 else:
#                     x=np.concatenate((x,xch),axis=0)
#                     y=np.concatenate((y,ych))
#                     z=np.concatenate((z,zch))
#     
#     
#     return x, y, z
# 
# 
# 
# 
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
# def get_child_particles(x_part,y_part,z_part,clumpid,child):
#     xtemp=np.zeros(len(x_part))
#     ytemp=np.zeros(len(x_part))
#     ztemp=np.zeros(len(x_part))
#     counter=0
#     for j in range(0,len(x_part)):
#         if (clumpid[j]==child):
#             xtemp[counter]=x_part[j]
#             ytemp[counter]=y_part[j]
#             ztemp[counter]=z_part[j]
#             counter+=1
# 
#     x=np.zeros(counter)
#     y=np.zeros(counter)
#     z=np.zeros(counter)
#     for j in range(counter):
#         x[j]=xtemp[j]
#         y[j]=ytemp[j]
#         z[j]=ztemp[j]
#     
#     ##########################################
#     # print out particles if needed
#     #         for j in range(len(x)):
#     #             print ('{0:12.7f}{1:12.7f}{2:12.7f}{3:12d}'.format(x[j], y[j], z[j], i+1))
#     ##########################################
# 
#     return x, y, z
# 
# 
# 
# 
# 
# 
# #----------------------------------------------------------

