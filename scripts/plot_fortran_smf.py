#!/usr/bin/env python3

#!/apps/vendor/Anaconda/4.3.1/anaconda3/bin/python3
#--------------------------------------------------
# Plots the stellar mass functions obtained from
# get_smf.f03. Run in directory, where
# output_XXXXX directories are in.
#--------------------------------------------------


import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                                'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)

# for legend
from matplotlib.font_manager import FontProperties # for legend
fontP=FontProperties()
fontP.set_size('x-small')
fontP.set_family('serif')




nfiles = 0
zlow = 0
zhigh = 0
zmean = 0
shortname = 0
obsdatalist = []
outputlist = []
#  obsdir = '/home/cluster/mivkov/masterarbeit/files/observational_data/behroozi-2013-data-compilation/smf_ms/'
obsdir = '/home/mivkov/UZH/Masterarbeit/masterarbeit/files/observational_data/behroozi-2013-data-compilation/smf_ms/'

# how many plots?
rows_smf = 6
cols_smf = 4

# snapshot data
redshift = 0
unit_l = 0
masses = 0
datalist = []
iz0 = 0

# physics
Mpc = 3.086e24 # Mpc/cm







#============================
def get_snapshot_data():
#============================
    """
    Reads in smf.txt and info_XXXXX.txt files.
    """

    from os import getcwd
    from os import listdir

    global redshift, unit_l, datalist, masses, outputlist

    workdir = getcwd()
    filelist = listdir(workdir)

    outputlist = []
    outdirtemp = 'output_'
    for filename in filelist:
        if filename.startswith(outdirtemp):
            # first check whether directory contains galaxy files
            txtfilelist = listdir(filename)
            #  if 'smf-new.txt' in txtfilelist:
            if 'smf-new.txt' in txtfilelist: # post 30.04.2021
                outputlist.append(filename)

    if len(outputlist)<1:
        print("I didn't find any output_XXXXX directories in current working directory.")
        print("Or none that contain an output_XXXXX/galaxies_XXXXX.txt00001 file.")
        print("Are you in the correct workdir?")
        quit()

    outputlist.sort()
    noutput = len(outputlist)
    print("First output containing galaxies: ", outputlist[0])

    a_exp = np.zeros(noutput, dtype='float')
    unit_l = np.zeros(noutput, dtype='float')

    #-----------------------------
    # Read info files
    #-----------------------------

    for i,out in enumerate(outputlist):
        fn = out+'/info_'+out[-5:]+'.txt'
        try:
            infofile = open(fn)
            for j in range(9):
                infofile.readline() # skip first 9 lines

            # get expansion factor
            aline = infofile.readline()
            astring, equal, aval = aline.partition("=")
            afloat = float(aval)
            a_exp[i] = afloat

            for j in range(5):
                infofile.readline() # skip first 9 lines

            lline = infofile.readline()
            lstring, equal, lval=lline.partition('=')
            lfloat = float(lval)
            unit_l[i] = lfloat/Mpc


        except IOError: # If file doesn't exist
            print("Didn't find any info data in ", srcdir)
            break

    redshift = 1.0/a_exp - 1




    #-------------------------------
    # Read smf.txt files
    #-------------------------------

    datalist = [0 for i in outputlist]
    for i,out in enumerate(outputlist):
        #  fn = out+'/smf.txt'
        fn = out+'/smf-new.txt'
        datalist[i] = np.loadtxt(fn, dtype='float', skiprows=1, usecols=([1]))

    # read masses only once
    #  masses = np.loadtxt(outputlist[0]+'/smf.txt', dtype='float', skiprows=1, usecols=[0])
    masses = np.loadtxt(outputlist[0]+'/smf-new.txt', dtype='float', skiprows=1, usecols=[0])

    return





#======================
def get_obs_data():
#======================
    """
    Read in observational data.
    """

    global nfiles, zlow, zhigh, zmean, shortname, obsdatalist

    zlow, zhigh, zmean, is_linear, obsfile, shortname = np.genfromtxt(obsdir+'contents.txt',
        comments="#", unpack=True, dtype=np.str)

    zlow = zlow.astype(np.float)
    zhigh = zhigh.astype(np.float)
    zmean = zmean.astype(np.float)
    zmean = zmean.astype(np.float)
    is_linear = is_linear.astype(np.int)

    nfiles = zlow.shape[0]

    obsdatalist = [0 for i in range(nfiles)]

    for i,f in enumerate(obsfile):
        obsdatalist[i] = np.loadtxt(obsdir+f, comments="#", dtype='float')

    # transform linear data to logarithmic data
    for i,l in enumerate(is_linear):
        if l > 0:
            d = obsdatalist[i]
            for f in range(d.shape[0]):
                temp =  np.log10(d[f,1])
                d[f,2] = np.log10(d[f,1]+d[f,2]) -temp
                if d[f,1]-d[f,3] < 0:
                    print("HEY", obsfile[i], d[f,1], d[f,3])
                d[f,3] = temp-np.log10(d[f,1]-d[f,3])
                d[f,1] = temp

    return




#===============================
def plot_smf():
#===============================
    """
    Create figure, plot obtained and observational SMFs.
    """


    #---------------------------------------------
    # find unique plots of observational data
    #---------------------------------------------

    uniques = 0
    zlow_plot  = np.ones(nfiles)
    zhigh_plot = np.ones(nfiles)
    zmean_plot = np.ones(nfiles)
    ax_index = np.ones(nfiles, dtype=int)
    zlow_plot  *= -1
    zhigh_plot *= -1
    zmean_plot *= -1


    for f in range(nfiles):
        for i in range(uniques+1):
            if (zlow_plot[i]==-1):
                zlow_plot[i] = zlow[f]
                zhigh_plot[i] = zhigh[f]
                zmean_plot[i] = zmean[f]
                ax_index[f] = i
                uniques += 1
                break
            elif ((zlow_plot[i]==zlow[f]) and (zhigh_plot[i]==zhigh[f])):
                ax_index[f] = i
                break




    #--------------------
    # Create figure
    #--------------------

    fig = plt.figure(figsize=(12,16),dpi=300)
    #  fig = plt.figure(figsize=(20,20))
    axes = [0 for i in range(rows_smf*cols_smf)]

    for i in range(uniques):
        axes[i] = fig.add_subplot(rows_smf, cols_smf, i+1)

    clist = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']



    #----------------------------
    # Plot observational data
    #----------------------------

    titles = ["" for f in range(uniques)]
    nlines = [0 for f in range(uniques)]

    for f in range(nfiles):
        i = ax_index[f]
        ax = axes[i]
        ax.errorbar(obsdatalist[f][:,0], obsdatalist[f][:,1],
            yerr=(obsdatalist[f][:,2], obsdatalist[f][:,3]),
            label=shortname[f], c=clist[2+nlines[i]]
            )
        nlines[i] += 1

        title=str(zlow_plot[ax_index[f]])+r'$< z <$'+str(zhigh_plot[ax_index[f]])
        titles[ax_index[f]] = title





    #----------------------------
    # Plot my data
    #----------------------------

    masses_plot = 0.5*(masses[1:]+masses[:-1])
    dex_divide = masses[1:]-masses[:-1]

    nbins=masses.shape[0]-1
    sm_main = np.zeros(nbins, dtype=float)

    # Note: smf are conventionally defined over comoving volume
    vol = (unit_l[0]*(1+redshift[0]))**3





    #-----------------------------------
    # 1) Plot smf over redshift interval
    #-----------------------------------

    for i in range(uniques):
        zh = zhigh_plot[i]
        zl = zlow_plot[i]
        sm_main[:] = 0
        volcount = 0

        for j,z in enumerate(redshift):
            if (z>=zl) and (z<=zh):
                sm_main[:] += datalist[j][1:][:]
                volcount += 1

        if (volcount>0):

            mask = sm_main>0
            dd = dex_divide[mask]
            sm = sm_main[mask]
            sm = np.log10(sm/vol/volcount/dd)

            if volcount > 1:
                stext = ' samples'
            else:
                stext = ' sample'

            axes[i].plot(masses_plot[mask], sm, c=clist[0],
                label='averaged $\Phi$ over z interval, '+str(volcount)+stext)




    #-------------------------------------------
    # 2) Plot smf at center of redshift interval
    #-------------------------------------------

    for i in range(uniques):
        zm = zmean_plot[i]
        sm_main[:] = 0

        j = np.argmin(np.absolute(redshift-zm)) # find output closest to z_mean_plot
        sm_main[:] = datalist[j][1:][:]

        mask = sm_main>0
        dd = dex_divide[mask]
        sm = sm_main[mask]
        sm = np.log10(sm/vol/dd)
        axes[i].plot(masses_plot[mask], sm, "--", c=clist[1],
            label='$\Phi$ at $z={0:5.3f}$'.format(abs(redshift[j])) )
        axes[i].legend()




    #---------------------------
    # Tweak plots
    #---------------------------

    from matplotlib.font_manager import FontProperties # for legend

    fontP=FontProperties()
    fontP.set_size('x-small') # sizes = ['xx-small', 'x-small', 'small', 'medium', 'large','x-large', 'xx-large']

    for r in range(rows_smf):
        for c in range(cols_smf):

            ind = c + cols_smf*r

            if ind<uniques:

                ax = axes[ind]

                ax.legend(title=titles[ind], loc='lower left', prop=fontP)
                ax.set_ylim(-6.2,-0.2)
                ax.set_xlim(5,13)
                #  ax.grid()
                ax.tick_params(axis='both', direction='in')
                ax.tick_params(axis='both', bottom=True, top=True, left=True, right=True)

                if c > 0:
                    ax.set_yticklabels([])
                if r > 0:
                    ax.set_xticklabels([])
                if r == 0:
                    ax.xaxis.tick_top()


    plt.subplots_adjust(left=0.06, right=0.99, top=0.96, bottom=0.02,wspace=0., hspace=0.)
    plt.figtext(0.02, 0.5, r"$\log_{10} \Phi(M_*)$ $[($Mpc$(1+z))^{-3}$ dex$^{-1}$]", ha='center', va='center',rotation=90, size=16)
    plt.figtext(0.5, 0.99, r"$\log_{10}$ $ M_* / M_{\odot}$ ", ha='center', va='center', size=16)

    #  plt.savefig('smf.pdf', format='pdf')
    plt.savefig('smf.png', dpi=300)


#==============
def main():
#==============

    get_obs_data()
    get_snapshot_data()

    plot_smf()








if __name__ == "__main__":
    main()