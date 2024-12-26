#!/usr/bin/python3
# test on which processor IDs particles are by finding particle id

import numpy as np
from sys import argv #command line arguments

#============================================================
def read_particles(start_ind, noutput):
#============================================================
    """
    extracts all particle data from files. 
    start_ind: index from which to start out from when reading filenames given as cmd line args
    particles:  number of particles of sim in total
    noutput:        number of files to go through
    
    """


    print("Reading in particle data.")
    

    for i in range(0,noutput):
        inputfile=str(argv[i+start_ind]) 
        temp_data=np.loadtxt(inputfile, dtype='float', skiprows=1, usecols=[0,1,2,6,8])
    
        counter = 0
        for p in temp_data[:,4]:
            if p<240000 and p>=200000:
                counter += 1

        print("Found ",counter, "particles on processor", i)

    return




read_particles(1, 4)
