#!/usr/bin/python3
#!/apps/vendor/Anaconda/4.3.1/anaconda3/bin/python3
#---------------------------------------------------
# Plots results obtained from get_smf.f03 for the
# report. This script looks for the following
# directories:
#   galaxy-512
#   galaxy-512-100MPC
#---------------------------------------------------



import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

# use LaTeX text
from matplotlib import rc
#  rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
# for Palatino and other serif fonts use:
#  rc('font',**{'family':'serif','serif':['Computer Modern']})
#  rc('text', usetex=True)

rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                                'monospace': ['Computer Modern Typewriter'],
                                'size':16})
rc('text', usetex=True)

# for legend
from matplotlib.font_manager import FontProperties # for legend
fontP=FontProperties()
#  fontP.set_size('x-small')
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
rows_smf = 2
cols_smf = 3

# snapshot data
redshift = 0
unit_l = 0
masses = 0
datalist = []
iz0 = 0

# physics
Mpc = 3.086e24 # Mpc/cm

# figure
figs = []
axes = []
clist = []
ax_index = []

# plot parameters from observational data
zlow_plot = 0
zhigh_plot = 0
zmean_plot = 0
titles = 0






#==============================
def get_snapshot_data(srcdir):
#==============================
    """
    Reads in smf.txt and info_XXXXX.txt files.
    srcdir: directory string that containts output_XXXXX dirs
    """

    from os import listdir

    global redshift, unit_l, datalist, masses, outputlist

    filelist = listdir(srcdir)

    outputlist = []
    for filename in filelist:
        if filename.startswith('output_'):
            # first check whether directory contains galaxy files
            fullfname = srcdir+'/'+filename
            txtfilelist = listdir(fullfname)
            if 'smf.txt' in txtfilelist:
                outputlist.append(fullfname)

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
            print("Didn't find any info data in ", out)
            break

    redshift = 1.0/a_exp - 1




    #-------------------------------
    # Read smf.txt files
    #-------------------------------

    datalist = [0 for i in outputlist]
    for i,out in enumerate(outputlist):
        fn = out+'/smf.txt'
        datalist[i] = np.loadtxt(fn, dtype='float', skiprows=1, usecols=([1]))

    masses = np.loadtxt(outputlist[0]+'/smf.txt', dtype='float', skiprows=1, usecols=[0])

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
def plot_smf(label, lineindex):
#===============================
    """
    plot obtained SMFs.
    label       : string for plot label
    lineindex   : index of line for colour
    """

    #----------------------------
    # Plot my data
    #----------------------------

    masses_plot = 0.5*(masses[1:]+masses[:-1])
    dex_divide = masses[1:]-masses[:-1]

    nbins=masses.shape[0]-1
    sm_main = np.zeros(nbins, dtype=float)

    # Note: smf are conventionally defined over comoving volume
    vol = (unit_l[0]*(1+redshift[0]))**3


    #-------------------------------------------
    # Plot smf at center of redshift interval
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
        axes[i].plot(masses_plot[mask], sm, c=clist[lineindex],
            label=label)
        axes[i].legend(prop=fontP)

    return







#========================
def setup_figure():
#========================
    """
    Setup figure, plot observational data
    """

    global fig, axes, uniques, clist
    global zlow_plot, zhigh_plot, zmean_plot, ax_index, titles

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

    global figs
    figs = [0 for i in range(4)]
    for i in range(4):
        figs[i] = plt.figure(figsize=(16,8),dpi=300)
    #  fig = plt.figure(figsize=(20,20))
    axes = [0 for i in range(rows_smf*cols_smf*4)]

    fnr = 0
    for i in range(uniques):
        if i - fnr*(rows_smf*cols_smf) == rows_smf*cols_smf:
            fnr += 1

        axes[i] = figs[fnr].add_subplot(rows_smf, cols_smf, i+1-fnr*(rows_smf*cols_smf))

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
            label=shortname[f], c=clist[2+nlines[i]],
            zorder=0, capsize=2, elinewidth=1
            )
        nlines[i] += 1

        #  title=str(zlow_plot[ax_index[f]])+r'$< z <$'+str(zhigh_plot[ax_index[f]])
        title=r'$z \sim {0:4.2f}$'.format(zmean_plot[ax_index[f]])
        titles[ax_index[f]] = title
    return




#========================
def tweak_plots():
#========================
    """
    Tweak plots and save figure.
    """

    #---------------------------
    # Tweak plots
    #---------------------------

    for i,f in enumerate(figs):
        plt.figure(f.number) # set current figure
        for r in range(rows_smf):
            for c in range(cols_smf):

                ind = c + cols_smf*r + i * (rows_smf*cols_smf)

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


        plt.subplots_adjust(left=0.06, right=0.99, top=0.92, bottom=0.02,wspace=0., hspace=0.)
        plt.figtext(0.02, 0.5, r"$\log_{10} \Phi(M_*)$ $[($Mpc$(1+z))^{-3}$ dex$^{-1}$]", ha='center', va='center',rotation=90, size=22)
        plt.figtext(0.5, 0.98, r"$\log_{10}$ $ M_* / M_{\odot}$ ", ha='center', va='center', size=22)


        plt.savefig('smf_presentation_both_sims-'+str(i+1)+'.pdf', format='pdf')

    return







#==============
def main():
#==============

    get_obs_data()
    setup_figure()
    for i,case in enumerate([['galaxy-512-69MPC', r'\texttt{G69}'], ['galaxy-512-100MPC', r'\texttt{G100}']]):
        get_snapshot_data(case[0])
        plot_smf(case[1], i)
    tweak_plots()

    return









if __name__ == "__main__":
    main()
