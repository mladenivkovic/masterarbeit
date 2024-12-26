#!/usr/bin/python2


#===========================================================================================
# Count how many halos at z = 0 you have with mass > mass_threshold
#
#
# -------------
#   Usage:
# -------------
#   eval_tree.py 
#
#
#===========================================================================================




#=======================
# declare globals
#=======================
outputnrs = 0
a_exp = 0
mass = 0
npart = 0

#todo: temp
root_current = 0








#======================
class global_params:
#======================
    """
    An object to store all global parameters in, so you can only
    pass 1 object to functions that need it.

    """


    #======================
    def __init__(self):
    #======================
        """
        Initialises object. 
        """
        self.lastdir = ''           # last output directory
        self.noutput = 0            # number of outputs to work with
        self.ncpu = 0               # how many processors were used for simulation run

        self.workdir = ''           # current working directory
        self.prefix=''              # filename prefix

        self.verbose = False        # whether to print more details of what this program is doing

        # dictionnary of accepted keyword command line arguments
        self.accepted_args = {
            '-v' : self.set_verbose,
            '--verbose' : self.set_verbose
            } 


        return






    #=============================
    def read_cmdlineargs(self):
    #=============================
        """
        Reads in the command line arguments and stores them in the
        global_params object.
        """
        from sys import argv 

        nargs = len(argv)
        i = 1 # first cmdlinearg is filename of this file, so skip it

        while i < nargs:
            arg  = argv[i]
            arg = arg.strip()
            if arg in self.accepted_args.keys():
                self.accepted_args[arg]()
            else:
                print "I didn't recognize the argument '", arg, "'"
                quit()

            i+= 1

        return





    #==========================
    def get_output_info(self):
    #==========================
        """
        Read in the output info based on the files in the current
        working directory.
        Reads in last directory, ncpu, noutputs. 
        """

        from os import getcwd
        from os import listdir

        self.workdir = getcwd()
        filelist = listdir(self.workdir)
        
        outputlist = []
        for filename in filelist:
            if "output_" in filename:
                outputlist.append(filename)

        
        if len(outputlist)<1:
            print "I didn't find any output_XXXXX directories in current working directory."
            print "Are you in the correct workdir?"
            quit()

        outputlist.sort()
        self.lastdir = outputlist[-1]
        self.lastdirnr=int(self.lastdir[-5:])
        self.noutput = len(outputlist)

        # read ncpu from infofile in last output directory
        infofile = self.lastdir+'/'+'info_'+self.lastdir[-5:]+'.txt'
        f = open(infofile, 'r')
        ncpuline = f.readline()
        line = ncpuline.split()
        
        self.ncpu = int(line[-1])
        f.close()

        return


    



    #========================
    # Setter methods
    #========================

    def set_verbose(self):
        self.verbose = True
        return








#============================================
def read_mergertree_data(params):
#============================================
    """
    reads in mergertree data as written by the mergertree patch.
    NOTE: requires mergertree.txtYYYYY files.

    parameters:
        params:     class global_params object, containing the global parameters

    returns:
        mass:      numpy array containing descendant masses
        npart:     nump array containing descendant particle numbers

    """ 

    import numpy as np
    import warnings

    noutput = params.noutput
    ncpu = params.ncpu

    if params.verbose:
        print "Reading in mergertree data."

    
    # create lists where to store stuff
    fname = 'mergertree.txt'


    startnr=params.lastdirnr
    outputnrs = np.array(range(startnr, startnr-noutput, -1))
    a_exp = np.zeros(noutput)
    unit_l = np.zeros(noutput)
    unit_m = np.zeros(noutput)
    unit_t = np.zeros(noutput)
    unit_dens = np.zeros(noutput)

    dir_template = 'output_'


    for output in range(noutput):

        dirnr =  str(startnr - output).zfill(5)
        srcdir = dir_template + dirnr

        try:
            #------------------------------------------------------
            # get time, redshift, and units even for output_00001
            #------------------------------------------------------
            fileloc = srcdir+'/info_'+dirnr+'.txt'
            infofile = open(fileloc)
            for i in range(9):
                infofile.readline() # skip first 8 lines
                # 9 if you don't read in time
            
            # get expansion factor
            aline = infofile.readline()
            astring, equal, aval = aline.partition("=")
            afloat = float(aval)
            a_exp[output] = afloat

            for i in range(5):
                infofile.readline() # skip 5 lines

            # get unit_l
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            unit_l[output] = unitfloat
            
            # get unit_dens
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            unit_dens[output] = unitfloat

            # get unit_t
            unitline = infofile.readline()
            unitstring, equal, unitval = unitline.partition("=")
            unitfloat = float(unitval)
            unit_t[output] = unitfloat
    
        except IOError: # If file doesn't exist
            print "Didn't find any info data in ", srcdir
            break


    
    redshift = 1/a_exp - 1 
    output_zero_ind = np.argmin(np.absolute(redshift))
    print "starting output is", outputnrs[output_zero_ind]





    dirnr =  str(outputnrs[output_zero_ind]).zfill(5)
    srcdir = dir_template + dirnr

    try:

        fnames = [srcdir + '/' + fname + str(cpu+1).zfill(5) for cpu in range(ncpu)]

        datalist = [np.zeros((1,3)) for i in range(ncpu)]
        i = 0
        for f in fnames:
            with warnings.catch_warnings():
                warnings.filterwarnings('error') # treat warnings as errors so I can catch them
                try:
                    datalist[i] = np.atleast_1d(np.loadtxt(f, dtype='float', skiprows=1, usecols=([3,4])))
                    i += 1
                except Warning:
                    continue

        if i > 0:
            fulldata = np.concatenate(datalist[:i], axis=0)
            mass = fulldata[:,0]
            npart = fulldata[:,1]


    except IOError: # If file doesn't exist
        print "Didn't find any progenitor data in", srcdir






    #--------------------------------------
    # Transform units to physical units
    #--------------------------------------

    Mpc = 3.086e24 # cm
    M_Sol = 1.98855e33 # g
    Gyear = (24*3600*365*1e9) # s

    unit_m = unit_dens*unit_l**3/M_Sol
    unit_l /= Mpc
    unit_t /= Gyear


    mass *= unit_m[output_zero_ind]



    return mass, npart









#====================
def main():
#====================

    import numpy as np

    # declare globals
    global outputnrs
    global mass
    global npart
    global mthresh



    #==================
    # Set up
    #==================
    params = global_params()
    params.read_cmdlineargs()
    params.get_output_info()


    

    #------------------
    # read in data
    #------------------
    mass, npart = read_mergertree_data(params)



    #---------------------------------------------
    # Set mass threshold for evaluations
    #---------------------------------------------

    mthresh = 5e11 # 10^12 M_Sol
    #  mthresh = 1e12 # 10^12 M_Sol
    # find how many particles that is

    for i in range(mass.shape[0]):
        if mass[i] > 0:
            mp = mass[i]/npart[i]
            mtstring = '{0:6E}'.format(mthresh)
            mpstring = '{0:6E}'.format(mp)
            print "mass threshold is "+mtstring+" solar masses, which corresponds to ~", int(mthresh/mp+0.5), "particles."
            print "particle mass is: "+mpstring+" solar masses" 
            break

    count = 0
    for m in mass:
        if m >= mthresh:
            count += 1

    print "Found", count, "clumps at z = 0 for mass threshold ", mthresh, "M_Sol"








    return






#=================================
if __name__ == "__main__":
#=================================
    main()
