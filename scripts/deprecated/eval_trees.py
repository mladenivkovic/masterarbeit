#!/usr/bin/python3


#===========================================================================================
# This script evaluates merger trees.
# Essentially gets histograms as is done in Sussing Merger Trees Project, and writes
# them to file so they can be plotted en masse with eval_plots.py
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

#  progenitors :           lists of numpy arrays of progenitors and descendants,
#  descendants :           starting with the last output step.
#  progenitor_outputnrs:   list of numpy arrays of the output number at which the
#                          progenitor is
#  outputnrs   :           numpy array the output number at which descendants were
#                          taken from
#  a_exp       :           numpy array of expansion factor correspondig to each
#                          output step
#  redshift    :           numpy array of redshifts corresponding to each output step
#  x, y, z     :           list of numpy arrays containing descendant positions
#  vx, vy, vz  :           list of numpy arrays containing descendant velocities
#  mass, npart :           lists of numpy arrays containing descendant masses and
#                          number of particles

lastdir = '0'
lastdirnr = ''
noutput = 0
ncpu = 0
workdir = ''
prefix = ''
output_start = 0
nhalosmax = 0


time = 0
rho_crit = 0
redshift = 0
a_exp = 0

outputnrs = 0
descendants = 0
progenitors = 0
progenitor_outputnrs = 0
mass = 0
npart = 0
x = 0
y = 0
z = 0
vx = 0
vy = 0
vz  = 0

mass_growth = 0
mass_flucts = 0
mffree = 0
displacements = 0
displfree = 0
mthresh = 0

counts = 0


unit_l = 0
ozi = 0


desc_ind_next = 0
desc_ind_current = 0
root_current = 0
root_next = 0
main_current = 0
main_next = 0



jumpers_data = 0
jumpers_ind = 0
jumpfree = 0
njumptot = 0

nkilled = 0
npartkilled = 0



#todo: temp
root_current = 0




import numpy as np









#===============================================
def calc_eval(o1, i1, o2, i2, o3, i3, debug):
#===============================================
    """
    o1 : output index of descendant
    i1 : array index of descendant
    o2 : output index of prog
    i2 : array index of prog
    o3 : output index of descendant's descendant
    i3 : array index of descendant's descendant
    """

    # declare globals
    global mass_growth, mass_flucts, mffree
    global displacements
    global displfree

    #todo: temp
    global root_current


    #------------------
    def r200(m, rc):
    #------------------
        """
        compute r200 for given mass m and rho_crit rc
        """
        R = (3 * m /(4 * np.pi * rc))**(1./3)
        return R

    # calc mass growth

    if mass[o1][i1] > mthresh and mass[o2][i2] > mthresh:

        # calc displacements
        u = vx[o1][i1] + vx[o2][i2]
        v = vy[o1][i1] + vy[o2][i2]
        w = vz[o1][i1] + vz[o2][i2]
        dt = times[o1] - times[o2]

        l = 0.5*(unit_l[o1]+unit_l[o2])

        dx = x[o1][i1] - x[o2][i2]
        if dx > 0.5 *l:
            dx = dx - l
        elif dx < -0.5*l:
            dx = l + dx
        dy = y[o1][i1] - y[o2][i2]
        if dy > 0.5*l:
            dy = dy - l
        elif dy < -0.5*l:
            dy = l + dy
        dz = z[o1][i1] - z[o2][i2]
        if dz > 0.5*l:
            dz = dz - l
        elif dz < -0.5*l:
            dz = l + dz

        if debug:
            dx0 = dx
            dy0 = dy
            dz0 = dz

        dx -= 0.5 * u * dt
        dy -= 0.5 * v * dt
        dz -= 0.5 * w * dt

        num = np.sqrt(dx**2 + dy**2 + dz**2)
        denum = 0.5*(r200(mass[o2][i2], rho_crit[o2]) + r200(mass[o1][i1],rho_crit[o1]) + np.sqrt(u**2+v**2+w**2)*dt)

        disp = num/denum



        displacements[displfree] = disp
        displfree += 1


        mass_growth[o1][i1] = 2.0/np.pi*np.arctan((mass[o1][i1]-mass[o2][i2])*(times[o1]+times[o2])/((mass[o1][i1]+mass[o2][i2])*(times[o1]-times[o2])))


        # calc mass fluctuation
        if mass_growth[o3][i3] < 2.0 and i3 >= 0: #don't allow default i3=-1
            mass_flucts[mffree] = (mass_growth[o3][i3] - mass_growth[o1][i1])*0.5
            mffree += 1





        if debug:
            pass
            #  msg = 'calculating {0:7d}->{1:7d} at output {2:3d}; got MG {3:8.3f}, displ {4:8.3f}'.format(descendants[o1][i1], progenitors[o1][i1], outputnrs[o1], mass_growth[o1][i1], displacements[displfree-1])
            #  more=True
            #  #  if abs(mass_growth[o1][i1])>10:
            #  #      msg += " !"
            #  #      more = True
            #  if displacements[displfree-1] > 10:
            #      msg += "   !!"
            #      more = True
            #      print( num, denum)
            #  if more:
            #      msg += ' N: {0:6.3f}, D: {1:6.3f}, dx: {2:6.3f} {3:6.3f} {4:6.3f}'.format(num, denum, dx0, dy0, dz0)
            #  print( msg)

    return






#===========================
def clean_up_jumpers():
#===========================
    """
    clean up jumpers.
    """

    global jumpers_data
    global jumpers_ind
    global jumpfree
    global njumptot

    # define new dtype to store jumpers
    jumptype = np.dtype([('root','i4'),('is main','bool'), ('desc out', 'i4'), ('desc ind', 'i4'),
        ('ddesc out', 'i4'), ('ddesc ind', 'i4')])

    # root:         index of root of tree
    # is main:      whether this is main branch
    # desc out:     output index for descendant
    # desc ind:     index of descenant in array
    # ddesc out:    output index for descendant's descendant
    # ddesc ind:    index of descendant's descendant in array

    jumpers_data = np.zeros(nhalosmax*10, dtype=jumptype)
    jumpers_ind = [np.ones(nhalosmax, dtype='int')*(-1) for i in range(noutput)]
    jumpfree = 0
    njumptot = 0

    # mark jumpers (at time of merging) as such
    for out in range(noutput):
        for i in range(progenitors[out].shape[0]):
            if progenitors[out][i] < 0:
                njumptot += 1
                sind = outputnrs[0]-progenitor_outputnrs[out][i]-1
                try:
                    # found a jumper: make descendant of original entry in past snapshots = 0
                    ind = get_jumper_ind(sind, -progenitors[out][i])
                    descendants[sind][ind] = 0
                except IndexError:
                    print( "Index Error @ marking jumpers", sind, outputnrs[sind], progenitor_outputnrs[out][i], descendants[out][i], progenitors[out][i])

    return






#===========================
def find_initial_roots():
#===========================
    global jumpers_data, jumpers_ind
    global jumpfree
    global root_current, root_next
    global main_current, main_next
    global desc_ind_current, desc_ind_next

    #---------------------------
    # determine initial roots
    #---------------------------

    # get root array: Store root of current branch
    root_current = np.ones(nhalosmax, dtype='int')*(-1)     # store index of the root
    root_next = np.ones(nhalosmax, dtype='int')*(-1)
    main_current = np.zeros(nhalosmax, dtype='bool')        # store whether this is main branch. zero = False
    main_next = np.zeros(nhalosmax, dtype='bool')
    desc_ind_current = np.ones(nhalosmax, dtype='int')*(-1) # store index of descendants two snapshots later
    desc_ind_next = np.ones(nhalosmax, dtype='int')*(-1)



    for i in range(descendants[ozi].shape[0]):
        if descendants[ozi][i] > 0:

            if progenitors[ozi][i] < 0:
                # if jumper: store root data of jumper
                oind = outputnrs[0] - progenitor_outputnrs[ozi][i] - 1
                iind = get_jumper_ind(oind, abs(progenitors[ozi][i]))

                jumpers_ind[oind][iind]=jumpfree
                jumpers_data[jumpfree]['root'] = i
                jumpers_data[jumpfree]['is main'] = True
                jumpers_data[jumpfree]['desc out'] = ozi
                jumpers_data[jumpfree]['desc ind'] = i
                jumpfree += 1
            elif progenitors[ozi][i] > 0:
                root_current[i] = i
                main_current[i] = True


        elif descendants[ozi][i] < 0:
            ind = -1
            try:
                for j, d in enumerate(descendants[ozi]):
                    if -d == descendants[ozi][i]:
                        ind = j
                        break
                if ind >= 0:
                    root_current[i] = ind
                    counts[ind]['nbranches'] += 1
                else:
                    print( "found no main for", descendants[ozi][i])
            except TypeError:
                print( "TypeError", type(descendants[ozi]))


    return






#==============================
def friedman(axp_min):
#==============================
    """
    Integrate friedman equation to get table of
    expansion factors and times.
    Gives a in units of H0.
    See ramses/utils/f90/utils.f90/subroutine friedman
    """

    def dadtau(axp_tau):
        dadtau = axp_tau**3 * (omega_m + omega_l*axp_tau**3 + omega_k * axp_tau)
        return np.sqrt(dadtau)

    def dadt(axp_t):
        dadt = 1/axp_t * (omega_m + omega_l*axp_t**3 + omega_k*axp_t)
        return np.sqrt(dadt)


    omega_m     =  0.272
    omega_l     =  0.728
    omega_k     =  0.0

    epsilon = 1e-4

    axp_tau = 1.0 # expansion factor
    axp_t = 1.0   # expansion factor
    tau = 0       # conformal time
    t = 0         # look-back time

    a_out = np.zeros(1000000, dtype='float')
    t_out = np.zeros(1000000, dtype='float')
    #  tau_out = np.zeros(1000000, dtype='float')
    H_out = np.zeros(1000000, dtype='float')

    i = 0
    while axp_tau >= axp_min or axp_t >= axp_min:
        dtau = epsilon * axp_tau / dadtau(axp_tau)
        axp_tau_pre = axp_tau - dadtau(axp_tau)*dtau/2
        axp_tau = axp_tau - dadtau(axp_tau_pre)*dtau
        tau = tau - dtau

        dt = epsilon * axp_t / dadt(axp_t)
        axp_t_pre = axp_t - dadt(axp_t)*dt/2
        axp_t = axp_t - dadt(axp_t_pre)*dt
        t = t - dt

        if (abs((t - t_out[i])/t) >= 1e-5):
            a_out[i] = axp_tau
            H_out[i] = dadtau(axp_tau)/axp_tau
            t_out[i] = t
            #  tau_out[i] = tau

            i += 1

    a_out[i] = axp_tau
    t_out[i] = t
    #  tau_out[i] = tau
    H_out[i] = dadtau(axp_tau)/axp_tau
    i+=1

    a_out = a_out[:i]
    t_out = t_out[:i]
    #  tau_out = tau_out[:i]
    H_out = H_out[:i]


    return a_out, t_out, H_out






#=============================
def get_cosmo_data():
#=============================
    """
    Compute times and Hubble parameter given the expansion factor.
    Directly return rho_crit at given a in
    M_Sol/Mpc^3 and cosmological time at that point.
    """

    global times, redshift, ozi, rho_crit


    H0 =  70.40
    times = np.zeros(a_exp.shape[0], dtype='float')
    H = np.zeros(a_exp.shape[0], dtype='float')

    a_out, t_out, H_out = friedman(a_exp.min())

    i = 1
    for j,a in enumerate(a_exp):
        while (a_out[i] > a) and (i < a_out.shape[0]-1):
            i += 1
        times[j] = t_out[i]*(a-a_out[i-1])/(a_out[i]-a_out[i-1]) + \
            t_out[i-1]*(a-a_out[i])/(a_out[i-1]-a_out[i])
        H[j] = H_out[i]*(a-a_out[i-1])/(a_out[i]-a_out[i-1]) + \
            H_out[i-1]*(a-a_out[i])/(a_out[i-1]-a_out[i])


    redshift = 1.0/a_exp - 1

    #----------------------
    # find where to start
    #----------------------
    # find redshift value closest to 0
    # (you need values of later times (times z<0) to check for jumpers correctly)
    ozi = np.argmin(np.absolute(redshift))
    print( "starting output is", outputnrs[ozi])


    #-----------------------
    # get physical units
    #-----------------------
    alpha = H0 * 1e5/3.08e24*(365*24*3600*1e9)

    times *= 1.0/alpha # get times in Gyrs: times were calculated in units of H_0

    # turn time around: t=0 at start to get cosmic time
    timestart = -times[-1]+times[ozi]
    times += timestart


    H *= alpha # get H in Gyrs^-1
    G = 4.492e-15 #Mpc^3/(M_sol Gyr^2)
    rho_crit = 3 * H**2 / (8 * np.pi *G)


    return






#=========================================
def get_jumper_ind(oind, prog_to_find):
#=========================================
    """
    get index where jumper progenitor was merged
    """

    for i, p in enumerate(progenitors[oind]):
        if p == prog_to_find:
            return i

    raise IndexError("Didn't find prog "+str(prog_to_find)+" in output "+str(outputnrs[oind]))






#==========================
def get_output_info():
#==========================
    """
    Read in the output info based on the files in the current
    working directory.
    Reads in last directory, ncpu, noutputs.
    """

    from os import getcwd
    from os import listdir

    global workdir, lastdir, lastdirnr, ncpu, noutput

    workdir = getcwd()
    filelist = listdir(workdir)

    outputlist = []
    for filename in filelist:
        if "output_" in filename:
            outputlist.append(filename)


    if len(outputlist)<1:
        print( "I didn't find any output_XXXXX directories in current working directory.")
        print( "Are you in the correct workdir?")
        quit()

    outputlist.sort()
    lastdir = outputlist[-1]
    lastdirnr=int(lastdir[-5:])
    noutput = len(outputlist)

    # read ncpu from infofile in last output directory
    infofile = lastdir+'/'+'info_'+lastdir[-5:]+'.txt'
    f = open(infofile, 'r')
    ncpuline = f.readline()
    line = ncpuline.split()
    ncpu = int(line[-1])
    f.close()

    return






#====================================
def histogram_and_write_data():
#====================================
    """
    Histogram Data and write to file
    """

    global npartkilled
    global njumptot

    outfname = 'eval_trees.txt'
    outfile = open(outfname, 'w')

    out = 'starting output is {0:6d}\n'.format(outputnrs[ozi])
    outfile.write(out)

    for i in range(descendants[ozi].shape[0]):
        if mass[ozi][i] > 0:
            mp = mass[ozi][i]/npart[ozi][i]
            break

    out = 'mass threshold is {0:12.4E} solar masses, which corresponds to ~ {1:9d} particles\n'.format(mthresh, int(mthresh/mp+0.5))
    outfile.write(out)

    npartkilled = npartkilled[:nkilled+1]
    nthresh = 100
    nthreshcount = 0
    for n in npartkilled:
        if n >= nthresh:
            nthreshcount += 1
    out = 'Removed {0:7d} clumps from tree, as no descendants were found. Max particles of kiled: {1:6.3f}; median: {2:6.3f}; above threshold of {3:5d} : {4:7d}\n'.format(nkilled, npartkilled.max(), np.median(npartkilled), nthresh, nthreshcount)
    outfile.write(out)

    out = 'Total number of jumpers in simulation: {0:10d}\n'.format(njumptot)
    outfile.write(out)



    # truncate junk from arrays
    descs = descendants[ozi][descendants[ozi]>0]
    progs = progenitors[ozi][descendants[ozi]>0]
    mbl = outputnrs[ozi]-counts[descendants[ozi]>0]['end of branch']
    nb = counts[descendants[ozi]>0]['nbranches']
    m = mass[ozi][descendants[ozi]>0]
    npt = npart[ozi][descendants[ozi]>0]


    out='number of haloes at z=0: {0:9d}\n'.format(descs.shape[0])
    outfile.write(out)
    out='nhalopart min: {0:9d}\n'.format(int(npt.min()+0.5))
    outfile.write(out)
    maxpart = int(npt.max()+0.5)
    out='nhalopart max: {0:9d}\n'.format(maxpart)
    outfile.write(out)
    out='nhalopart med: {0:9d}\n'.format(int(np.median(npt)+0.5))
    outfile.write(out)



    # make histograms
    bins = np.array([0, 101, 501, 1001, maxpart+1])
    bincounts, binedges = np.histogram(npt, bins=bins)
    hist_mbl, binedges = np.histogram(npt, bins=bins, weights=mbl)
    hist_nb, binedges = np.histogram(npt, bins=bins, weights=nb)

    average_mbl = hist_mbl*1.0 / bincounts
    average_nb = hist_nb*1.0 / bincounts


    for i in range(bins.shape[0]-1):
        out = 'bin {2:3d} particles: average main branch length: {0:12.3f} ; average number of branches: {1:12.3f} \n'.format(average_mbl[i], average_nb[i], i)
        outfile.write(out)



    binind = np.digitize(npt, bins=bins) # index of bin edge to the right

    #---------------------
    # write redshift
    #---------------------
    for i, z in enumerate(redshift[ozi:noutput-ozi+1]):
        out = '{0:12.4f}'.format(z)
        outfile.write(out)
        if i < redshift[ozi:noutput-ozi+1].shape[0]-1:
            outfile.write(",")
    outfile.write("\n")

    #---------------------
    # write particle bins
    #---------------------
    for i,b in enumerate(bins):
        out = '{0:6d}'.format(b)
        outfile.write(out)
        if i < bins.shape[0]-1:
            outfile.write(",")
    outfile.write("\n")


    #----------------------------------
    # write length of main branch
    #----------------------------------

    maxlen = mbl.max()
    while maxlen%5 != 0:
        maxlen += 1
    bins2 = np.array(range(1, maxlen+1))
    for i in range(1,bins.shape[0]):
        mblens, binedges = np.histogram(mbl[binind==i], bins=bins2)
        mblens = mblens * 1.0
        mblens /= mblens.sum() # normalise
        write_array(outfile, mblens, binedges[:-1]-1)
        # move binedges back: histogram sorts by [bin i, bin i+1), so len 2 will be in bin 3


    #----------------------------------
    # write number of branches
    #----------------------------------
    maxnb = nb.max()
    while maxnb%5 != 0:
        maxnb += 1
    bins2 = np.array(range(1, maxnb+1))
    for i in range(1,bins.shape[0]):
        nbs, binedges = np.histogram(nb[binind==i], bins=bins2)
        nbs = nbs * 1.0
        nbs /= nbs.sum() # normalise
        write_array(outfile, nbs, binedges[:-1]-1)



    #  logscale = np.logspace(-1, 0.0, num=50)
    #  bins2 = np.concatenate((-1.0*logscale[::-1], logscale))
    #  bins2 = np.linspace(-1.0, 1.0, 400)
    bins2=100

    #---------------------
    # write mass growth
    #---------------------
    mg = np.concatenate([mg_out[mg_out<2.0] for mg_out in mass_growth])
    mgplot, binedges = np.histogram(mg, bins=bins2)
    mgplot = mgplot*1.0
    mgplot /= mgplot.sum() # normalize
    write_array(outfile, mgplot, (binedges[1:]+binedges[:-1])/2.0)


    #----------------------------
    # write mass fluctuations
    #----------------------------
    mf = mass_flucts[mass_flucts < 2.0]
    mfplot, binedges = np.histogram(mf, bins=bins2)
    mfplot = mfplot*1.0
    mfplot /= mfplot.sum() # normalize
    write_array(outfile, mfplot, (binedges[1:]+binedges[:-1])/2.0)



    #----------------------------
    # write displacements
    #----------------------------
    disps = displacements[displacements != 0]
    dispsplot, binedges = np.histogram(disps, bins=bins2)
    dispsplot = dispsplot * 1.0
    dispsplot /= dispsplot.sum() # normalize
    write_array(outfile, dispsplot, (binedges[1:]+binedges[:-1])/2.0)


    print( "Finished evaluation. Results are written to", outfname)

    return






#============================================
def read_mergertree_data():
#============================================
    """
    reads in mergertree data as written by the mergertree patch.
    NOTE: requires mergertree.txtYYYYY files.
    """

    import warnings
    import gc


    print( "Reading in mergertree data.")


    global progenitors, descendants, progenitor_outputnrs, outputnrs
    global mass, npart
    global x, y, z, vx, vy, vz
    global a_exp, unit_l


    # create lists where to store stuff

    progenitors = [np.zeros(1) for i in range(noutput)]
    descendants = [np.zeros(1) for i in range(noutput)]
    progenitor_outputnrs = [np.zeros(1) for i in range(noutput)]
    mass = [np.zeros(1) for i in range(noutput)]
    npart = [np.zeros(1) for i in range(noutput)]
    x = [np.zeros(1) for i in range(noutput)]
    y = [np.zeros(1) for i in range(noutput)]
    z = [np.zeros(1) for i in range(noutput)]
    vx = [np.zeros(1) for i in range(noutput)]
    vy = [np.zeros(1) for i in range(noutput)]
    vz = [np.zeros(1) for i in range(noutput)]
    outcount = 0

    startnr=lastdirnr
    outputnrs = np.array(range(startnr, startnr-noutput, -1))
    a_exp = np.zeros(noutput)
    unit_l = np.zeros(noutput)
    unit_m = np.zeros(noutput)
    unit_t = np.zeros(noutput)
    unit_dens = np.zeros(noutput)

    dir_template = 'output_'

    # define new datatype for mergertree output
    mtree = np.dtype([('clump', 'i4'), ('prog', 'i4'), ('prog_outnr', 'i4'),
        ('mass', 'f8'), ('npart', 'f8'),
        ('x', 'f8'), ('y', 'f8'), ('z', 'f8'),
        ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8')])


    #---------------------------
    # Loop over directories
    #---------------------------

    for output in range(noutput):
        # loop through every output: Progenitor data only starts at output_00002,
        # but you'll need time/redshift data from output_00001!

        # Start with last directory (e.g. output_00060),
        # work your way to first directory (e.g. output_00001)
        dirnr =  str(startnr - output).zfill(5)
        srcdir = dir_template + dirnr

        if output < noutput -1: # don't try to read progenitor stuff from output_00001

            # Stop early if you reach a directory that has no mergertree.txt* files
            # (Can happen if there are no halos in the simulation yet)
            try:
                warnings.simplefilter("ignore", ResourceWarning)
                warnings.simplefilter("ignore", UserWarning)


                #------------------------------
                # Read in mergertree data
                #------------------------------

                fnames = [srcdir + '/' + "mergertree_"+dirnr+'.txt' + str(cpu+1).zfill(5) for cpu in range(ncpu)]

                datalist = [np.zeros((1,3)) for i in range(ncpu)]
                i = 0
                for f in fnames:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('error') # treat warnings as errors so I can catch them
                        try:
                            datalist[i] = np.atleast_1d(np.genfromtxt(f, dtype=mtree, skip_header=1))
                            i += 1
                        except Warning:
                            continue


                #---------------------------------
                # Sort out data
                #---------------------------------

                if i > 0:
                    fulldata = np.concatenate(datalist[:i], axis=0)

                    descendants[output] = fulldata[:]['clump']
                    progenitors[output] = fulldata[:]['prog']
                    progenitor_outputnrs[output] = fulldata[:]['prog_outnr']
                    mass[output] = fulldata[:]['mass']
                    npart[output] = fulldata[:]['npart']
                    x[output] = fulldata[:]['x']
                    y[output] = fulldata[:]['y']
                    z[output] = fulldata[:]['z']
                    vx[output] = fulldata[:]['vx']
                    vy[output] = fulldata[:]['vy']
                    vz[output] = fulldata[:]['vz']

                outcount += 1

            except IOError: # If file doesn't exist
                print( "Didn't find any mergertree data in", srcdir)





        try:
            #------------------------------------------------------
            # get time, redshift, and units even for output_00001
            #------------------------------------------------------
            fileloc = srcdir+'/info_'+dirnr+'.txt'
            infofile = open(fileloc)
            for i in range(9):
                infofile.readline() # skip first 9 lines

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

            infofile.close()

        except IOError: # If file doesn't exist
            print( "Didn't find any info data in ", srcdir)
            break


    # keep only entries that contain data
    if outcount > 1:
        descendants = descendants[:outcount]
        progenitors = progenitors[:outcount]
        progenitor_outputnrs = progenitor_outputnrs[:outcount]
        mass = mass[:outcount]
        npart = npart[:outcount]
        x = x[:outcount]
        y = y[:outcount]
        z = z[:outcount]
        vx = vx[:outcount]
        vy = vy[:outcount]
        vz = vz[:outcount]




    #--------------------------------------
    # Transform units to physical units
    #--------------------------------------

    Mpc = 3.086e24 # cm
    M_Sol = 1.98855e33 # g
    Gyear = (24*3600*365*1e9) # s

    unit_m = unit_dens*unit_l**3/M_Sol
    unit_l /= Mpc
    unit_t /= Gyear

    # transform units to physical units
    for i in range(len(descendants)):
        x[i] *= unit_l[i] # only transform later when needed; Need to check for periodicity first!
        y[i] *= unit_l[i]
        z[i] *= unit_l[i]
        mass[i] *= unit_m[i]
        vx[i] *= unit_l[i]/unit_t[i]
        vy[i] *= unit_l[i]/unit_t[i]
        vz[i] *= unit_l[i]/unit_t[i]


    # collect garbage
    del unit_m, unit_t, unit_dens
    gc.collect()



    return







#======================
def walk_trees():
#======================
    """
    Create Trees and count/calculate stuff
    """

    global jumpers_data, jumpers_ind
    global jumpfree
    global root_current, root_next
    global main_current, main_next
    global desc_ind_current, desc_ind_next
    global nkilled, npartkilled

    # counters
    nkilled = 0
    npartkilled = np.zeros(nhalosmax*noutput)

    # debug
    root_debug = 114

    for out in range(ozi+1, noutput):
        for i in range(descendants[out].shape[0]):



            #----------------------
            # IF MAIN PROGENITOR
            #----------------------
            if descendants[out][i] > 0:

                # get index in later snapshot
                parent_ind = np.where(progenitors[out-1]==descendants[out][i])
                # check that you found something
                if len(parent_ind[0])<1:
                    nkilled += 1
                    npartkilled[nkilled]=npart[out][i]
                    continue
                else:
                    parent_ind = parent_ind[0][0]
                    desc_ind_next[i] = parent_ind


                if descendants[out-1][parent_ind] != 0:
                    calc_eval(out-1, parent_ind, out, i,
                        out-2, desc_ind_current[parent_ind],
                        root_current[parent_ind]==root_debug)
                else:
                    jind = jumpers_ind[out-1][parent_ind]
                    calc_eval(out-1, parent_ind, out, i,
                    jumpers_data[jind]['desc out'], jumpers_data[jind]['desc ind'],
                    False)

                # if this clump is in the tree at z=0:
                if root_current[parent_ind] >= 0:
                    # main progenitor found.
                    # if is in main branch, add count. Otherwise, don't count for length of side branches
                    if main_current[parent_ind]:
                        counts[root_current[parent_ind]]['end of branch'] = progenitor_outputnrs[out][i]


                    if progenitors[out][i] > 0:
                        # inherit root index and whether this is main branch
                        main_next[i] = main_current[parent_ind]
                        root_next[i] = root_current[parent_ind]

                    elif progenitors[out][i] < 0:
                        # if jumper: store root data of jumper
                        oind = outputnrs[0] - progenitor_outputnrs[out][i]-1
                        # write where progenitor = progenitor; descendant will be = 0, so won't be found
                        iind = get_jumper_ind(oind, -progenitors[out][i])
                        jumpers_ind[oind][iind] = jumpfree
                        jumpers_data[jumpfree]['root'] = root_current[parent_ind]
                        jumpers_data[jumpfree]['is main'] = main_current[parent_ind]
                        jumpers_data[jumpfree]['desc out'] = out
                        jumpers_data[jumpfree]['desc ind'] = i
                        jumpers_data[jumpfree]['ddesc out'] = out - 1
                        jumpers_data[jumpfree]['ddesc ind'] = parent_ind
                        jumpfree += 1



            #----------------------
            # IF MERGER
            #----------------------
            elif descendants[out][i] < 0:

                # get index in later snapshot
                parent_ind = np.where(progenitors[out-1]==-descendants[out][i])
                # check that you found something
                if len(parent_ind[0])<1:
                    nkilled += 1
                    npartkilled[nkilled]=npart[out][i]
                    continue
                else:
                    parent_ind = parent_ind[0][0]


                if descendants[out-1][parent_ind] != 0:
                    calc_eval(out-1, parent_ind, out, i,
                        out-2, desc_ind_current[parent_ind],
                        root_current[parent_ind]==root_debug)
                else:
                    jind = jumpers_ind[out-1][parent_ind]
                    calc_eval(out-1, parent_ind, out, i,
                    jumpers_data[jind]['desc out'], jumpers_data[jind]['desc ind'],
                    False)

                # if this clump is in the tree at z=0:
                if root_current[parent_ind] >= 0:
                    # merger found. Count new branch
                    counts[root_current[parent_ind]]['nbranches'] += 1
                    root_next[i] = root_current[parent_ind]




            #-------------------------
            # IF JUMPER
            #-------------------------
            else: # descendant = 0
                jind = jumpers_ind[out][i]
                if jumpers_data[jind]['root'] >= 0:
                    #  if progenitors[out][i] > 0: # this is a given
                    if jumpers_data[jind]['is main']:
                        # add a count
                        counts[root_current[parent_ind]]['end of branch'] = progenitor_outputnrs[out][i]
                        # calc mass growth
                        oind = jumpers_data[jind]['desc out']   # output index for descendant
                        iind = jumpers_data[jind]['desc ind']   # index of descenant in array
                        ooind = jumpers_data[jind]['ddesc out'] # output index for descendant's descendant
                        iiind = jumpers_data[jind]['ddesc ind'] # index of descendant's descendant in array

                        calc_eval(oind, iind, out, i,
                            ooind, iiind,
                            root_current[parent_ind]==root_debug)

                main_next[i] = jumpers_data[jind]['is main']
                root_next[i] = jumpers_data[jind]['root']



        # overwrite root/main current for next round
        root_current = np.copy(root_next)
        main_current = np.copy(main_next)
        desc_ind_current = np.copy(desc_ind_next)
        # reset values for next round
        root_next[:] = -1
        main_next[:] = False
        desc_ind_next[:] = -1



    #  main_branch_len = outputnrs[ozi] - counts[:]['end of branch']
    #
    #  for i in range(descendants[ozi].shape[0]):
    #      if descendants[ozi][i] > 0:
    #          print( i, "halo", descendants[ozi][i], "has length of main branch",)
    #          print( main_branch_len[i], " and ", counts[i]['nbranches'], "branches")




    return







#====================================================
def write_array(outfile, counts, binedges):
#====================================================
    """
    Writes histogrammed arrays to file.
    """

    for i, c in enumerate(counts):
        out = '{0:15.7f} '.format(c)
        outfile.write(out)
        if i < counts.shape[0]-1:
            outfile.write(",")
    outfile.write("\n")


    for i,b in enumerate(binedges):
        out = '{0:15.7f} '.format(b)
        outfile.write(out)
        if i < binedges.shape[0]-1:
            outfile.write(",")
    outfile.write("\n")


    return






#====================
def main():
#====================

    global descendants
    global mass_growth, mass_flucts, mffree
    global displacements
    global displfree
    global mthresh
    global output_start, nhalosmax, noutput
    global desc_ind_next, desc_ind_current, counts

    # todo: temp
    global root_current



    #==================
    # Set up
    #==================
    get_output_info()
    read_mergertree_data()
    get_cosmo_data()




    #==================
    # Preparations
    #==================

    #---------------------------------------------
    # Set mass threshold for evaluations
    #---------------------------------------------

    mthresh = 5e11 # M_Sol
    #  mthresh = 1e13
    #  mthresh = 1e12 # M_Sol
    # find how many particles that is

    for i in range(descendants[ozi].shape[0]):
        if mass[ozi][i] > 0:
            mp = mass[ozi][i]/npart[ozi][i]
            mtstring = '{0:6E}'.format(mthresh)
            mpstring = '{0:6E}'.format(mp)
            print( "mass threshold is "+mtstring+" solar masses, which corresponds to ~", int(mthresh/mp+0.5), "particles.")
            print( "particle mass is: "+mpstring+" solar masses" )
            break


    #--------------------------------------------
    # get actual number of outputs we work with
    #--------------------------------------------
    noutput = len(descendants)
    output_start = outputnrs[ozi]

    #--------------------------------------------
    # get maximal nr of haloes
    #--------------------------------------------
    nhalos = [descendants[i].shape[0] for i in range(len(descendants))]
    nhalosmax = max(nhalos)

    #--------------------------------------------
    # define datatype for results
    #--------------------------------------------
    restype = np.dtype([('end of branch', 'i4'), ('nbranches', 'i4')])

    #--------------------------------------------
    # initialize count array
    #--------------------------------------------
    counts = np.empty(descendants[ozi].shape[0], dtype=restype)
    counts[:]['nbranches'] = 1 # include main branch as a branch!
    counts[:]['end of branch'] = outputnrs[ozi+1] # add 1, so minimal existance will be 1 snapshot


    #---------------------------------------
    # initialise other arrays you need
    #---------------------------------------

    # max will be 1, so initialise to 2 to see where value was changed
    mass_flucts = np.ones(nhalosmax*noutput)*2.0
    mffree = 0 # index for mass flucts
    mass_growth = [np.ones(nhalosmax, dtype='float')*2.0 for i in range(noutput)]

    displacements = np.zeros(nhalosmax*noutput)
    displfree = 0





    #==================================
    # GET STUFF DONE
    #==================================

    clean_up_jumpers()
    find_initial_roots()
    walk_trees()
    histogram_and_write_data()

    return









#=================================
if __name__ == "__main__":
#=================================
    main()
