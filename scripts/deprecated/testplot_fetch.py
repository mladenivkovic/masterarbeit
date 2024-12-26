# A module for the ramsesplot scripts
# Containts all functions about getting data

# Contains:
# get_halo_particles
# get_all_level_particles
# read_particles


import numpy as np
from sys import argv #command line arguments
import warnings


#########################
##--- ACQUIRING DATA   ##
#########################




#===================================
def find_halo(start_ind, noutput):
#===================================

    """
    Finds halo for a simulation assumed to have only 1.

    start_ind:  index where clump_* files start on cmd arg list
    noutput:    number of files through which to look through

    returns:
        halo:       halo ID of said halo.
    """

    # Get list of children for halo
    # noutput: number of outputfiles


    # read in clump data
    import ramsesplot_fetch as rf
    data_int, data = rf.get_clump_data(start_ind, noutput)
    
    
    #identify all children and their levels recursively for given halo
    i=0
    while(i<data_int.shape[0]):
        try:
            if (data_int[i,0] == data_int[i,2]): 
                halo = data_int[i,0]
                return halo 
            else:
                i+=1
        except IndexError:
            print("Caught exception!", data_int)
            halo = data_int[i]
            return halo


    print("Something went wrong in tesplot_fetch.find_halo(). Exiting.")
    quit()
    return



#===========================================
def get_all_clump_data(start_ind, noutput):
#===========================================

    """
    finds all halos and subhalos.

    start_ind:  index where clump_* files start on cmd arg list
    noutput:    number of files through which to look through

    returns:
        halos:          list of halo IDs.
        subhalos:       list of subhalo IDs.
    """

    # Get list of children for halo
    # noutput: number of outputfiles


    # read in clump data
    import ramsesplot_fetch as rf
    data_int, data = rf.get_clump_data(start_ind, noutput)
   

    halos = []
    subhalos = []
    
    #identify all children and their levels recursively for given halo

    try:
        for i in range(data_int.shape[0]):
            if (data_int[i,0] == data_int[i,2]): 
                halos.append(data_int[i,0])
            else:
                subhalos.append(data_int[i,0])
    except IndexError:
        # if data_int isn't a 2D array, there is only
        # 1 entry for the entire output.
        # That must be a halo
        print("Caught exception!", data_int)
        halos = [data_int[i]]

    

    return halos, subhalos







#======================================================================
def get_all_clump_particles(noutput, particles, halos, subhalos):
#======================================================================
    """
    Sorts all particles by level.

    noutput:        number of files to go through
    particles:      number of particles of sim in total
    subhalos:       list of all subhalos
    halos:          list of all halos


    returns:
        halox, haloy, haloz, subhalox, subhaloy, subhaloz
        lists of lists of x,y,z positions of particles that belong to halo/subhalo
        at same place of halos/subhalos lists
    """

    # read in particles
    data = read_particles(4+noutput, particles, noutput)

    print("Shape of read data:", data.shape)
    

    halox = []
    haloy = []
    haloz = []
    subhalox = []
    subhaloy = []
    subhaloz = []


    for halo in halos:
        halox.append([])
        haloy.append([])
        haloz.append([])

    for subhalo in subhalos:
        subhalox.append([])
        subhaloy.append([])
        subhaloz.append([])





    print("Sorting particles by halos.")

    for i in range(data.shape[0]):
        clumpID = int(abs(data[i,3]))
        if (clumpID in halos):
            haloind = halos.index(clumpID)
            halox[haloind].append([data[i,0]])
            haloy[haloind].append([data[i,1]])
            haloz[haloind].append([data[i,2]])

        elif (clumpID in subhalos): 
            subhaloind = subhalos.index(clumpID)
            subhalox[subhaloind].append([data[i,0]])
            subhaloy[subhaloind].append([data[i,1]])
            subhaloz[subhaloind].append([data[i,2]])
          
         
    return halox, haloy, haloz, subhalox, subhaloy, subhaloz




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
        px, py, pz, hx, hy, hz, mbx, mby, mbz, mbhx, mbhy, mbhz

    px, py, pz: list of numpy arrays of satellite halo particles, aranged by level
    hx, hy, hz: list of numpy arrays numpy array of halo particles, aranged by level
    mbx, mby, mbz: list of numpy arrays of most bound particles of satellite halos, aranged by level 
    mbhx, mbhy, mbhz : list of numpy arrays of most bound particles of main halos

    """

    # read in particles
    data = read_particles(3+noutput, particles, noutput)

    print("Shape of read data:", data.shape)
    

    



    # filter out all non-clump particles
    # store child particles in array for each child.
    # chose all particles for first guess of length.
    
    maxlevel = len(clump_levels)

    halox=np.zeros(particles)
    haloy=np.zeros(particles)
    haloz=np.zeros(particles)
    
    halombx = np.zeros(particles)
    halomby = np.zeros(particles)
    halombz = np.zeros(particles)

   

    parx = []
    pary = []
    parz = []

    parmbx = []
    parmby = []
    parmbz = []
    for l in range(maxlevel):
        parx.append(np.zeros(particles))
        pary.append(np.zeros(particles))
        parz.append(np.zeros(particles))
        parmbx.append(np.zeros(particles))
        parmby.append(np.zeros(particles))
        parmbz.append(np.zeros(particles))


    # set up found particle counter
    hc=0
    hmbc = 0
    pc=[0]*maxlevel
    pmbc=[0]*maxlevel




    print("Sorting particles by levels.")

    counter = 0
    for p in data[:,3]:
        if p == 0:
            counter += 1

    print("particles not in a clump:", counter)

    for i in range(data.shape[0]):
        clumpID = data[i,3]
        if (clumpID > 0): #particle is in clump
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

        elif (clumpID < 0): #particle is a most bound one           
            if -clumpID in halos:
                halombx[hmbc]=data[i,0]
                halomby[hmbc]=data[i,1]
                halombz[hmbc]=data[i,2]
                hmbc+=1
            else:
                for level in range(maxlevel):
                    if -clumpID in clump_levels[level]:
                        parmbx[level][pmbc[level]]=data[i,0]
                        parmby[level][pmbc[level]]=data[i,1]
                        parmbz[level][pmbc[level]]=data[i,2]
                        pmbc[level]+=1
                        break

             

    hx=np.zeros(hc)
    hy=np.zeros(hc)
    hz=np.zeros(hc)

    for i in range(hc):
        hx[i]=halox[i]
        hy[i]=haloy[i]
        hz[i]=haloz[i]

    mbhx=np.zeros(hmbc)
    mbhy=np.zeros(hmbc)
    mbhz=np.zeros(hmbc)

    for i in range(hmbc):
        mbhx[i]=halombx[i]
        mbhy[i]=halomby[i]
        mbhz[i]=halombz[i]
       
    
    px = []
    py = []
    pz = []

    mbx = []
    mby = []
    mbz = []

    
    for l in range(maxlevel):
        px.append(np.zeros(pc[l]))
        py.append(np.zeros(pc[l]))
        pz.append(np.zeros(pc[l]))
        mbx.append(np.zeros(pmbc[l]))
        mby.append(np.zeros(pmbc[l]))
        mbz.append(np.zeros(pmbc[l]))
        
        for p in range(pc[l]):
            px[l][p] = parx[l][p]
            py[l][p] = pary[l][p]
            pz[l][p] = parz[l][p]                   

        for p in range(pmbc[l]):
            mbx[l][p] = parmbx[l][p]
            mby[l][p] = parmby[l][p]
            mbz[l][p] = parmbz[l][p]
        
    return px, py, pz, hx, hy, hz, mbx, mby, mbz, mbhx, mbhy, mbhz




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
        fx, fy, fz, mux, muy, muz, hx, hy, hz, hmbx, hmby, hmbz

    fx, fy, fz:         lists of numpy arrays of particles. 
                        fx[i] is the particle list for child children[i]
    mbx, mby, mbz:      lists of most bound particles for child
    hx, hy, hz:         numpy array of halo particles.
    hmbx, hmby, hmbz:   lists of most bound particles for halo

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

    filteredmbx = []
    filteredmby = []
    filteredmbz = []
    

    for c in range(nchildren):
        filteredx.append(np.zeros(particles))
        filteredy.append(np.zeros(particles))
        filteredz.append(np.zeros(particles))

        filteredmbx.append(np.zeros(particles))
        filteredmby.append(np.zeros(particles))
        filteredmbz.append(np.zeros(particles))

    halox=np.zeros(particles)
    haloy=np.zeros(particles)
    haloz=np.zeros(particles)
    
    halombx = np.zeros(particles)
    halomby = np.zeros(particles)
    halombz = np.zeros(particles)

    # set up found particle counters
    fc=nchildren*[0]
    fmbc=nchildren*[0]
    hc=0
    hmbc=0

    for i in range(particles):
        if (len(children)>0):
        #find children particles
            if data[i,3] in children:
                childind = children.index(data[i,3])
                # if particle belongs to child, append to correct array
                filteredx[childind][fc[childind]]=data[i,0]  
                filteredy[childind][fc[childind]]=data[i,1]
                filteredz[childind][fc[childind]]=data[i,2] 
                # raise found particle index
                fc[childind]+=1

            elif -1*data[i,3] in children:
                childind = children.index(-1*data[i,3])
                # if particle belongs to child, append to correct array
                filteredmbx[childind][fmbc[childind]]=data[i,0]  
                filteredmby[childind][fmbc[childind]]=data[i,1]
                filteredmbz[childind][fmbc[childind]]=data[i,2] 
                # raise found particle index
                fmbc[childind]+=1


        if (data[i,3]==halo):
            #find halo particles
            # append halo particles 
            halox[hc]=data[i,0]
            haloy[hc]=data[i,1]
            haloz[hc]=data[i,2]
            hc+=1

        elif (-1*data[i,3]==halo):
            #find halo particles
            # append halo particles 
            halombx[hmbc]=data[i,0]
            halomby[hmbc]=data[i,1]
            halombz[hmbc]=data[i,2]
            hmbc+=1
           

   

    #create empty lists for final particle lists
    fx = []
    fy = []
    fz = []
    mbx = []
    mby = []
    mbz = []

    #add array for each child with correct number of particles
    for c in range(nchildren):
        fx.append(np.zeros(fc[c]))
        fy.append(np.zeros(fc[c]))
        fz.append(np.zeros(fc[c]))

        mbx.append(np.zeros(fmbc[c]))
        mby.append(np.zeros(fmbc[c]))
        mbz.append(np.zeros(fmbc[c]))

    #copy the particles into final arrays
    for c in range(nchildren):
        for i in range(fc[c]):
            fx[c][i]=filteredx[c][i]
            fy[c][i]=filteredy[c][i]
            fz[c][i]=filteredz[c][i]

        for i in range(fmbc[c]):
            mbx[c][i] = filteredmbx[c][i]
            mby[c][i] = filteredmby[c][i]
            mbz[c][i] = filteredmbz[c][i]



    hx=np.zeros(hc)
    hy=np.zeros(hc)
    hz=np.zeros(hc)

    hmbx=np.zeros(hmbc)
    hmby=np.zeros(hmbc)
    hmbz=np.zeros(hmbc)

    for i in range(hc):
        hx[i]=halox[i]
        hy[i]=haloy[i]
        hz[i]=haloz[i]

    for i in range(hmbc):
        hmbx[i] = halombx[i]
        hmby[i] = halomby[i]
        hmbz[i] = halombz[i]




    return fx, fy, fz, mbx, mby, mbz, hx, hy, hz, hmbx, hmby, hmbz




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



