#!/usr/bin/python3
#---------------------------------------------------
# Plots results obtained from get_smf.f03 for the
# report. This script looks for the following
# directories:
#   galaxy-512-new-cosmo-69MPC
#   galaxy-512-new-cosmo-100MPC
#---------------------------------------------------



import numpy as np
import matplotlib as mpl
#  mpl.use('Agg')
import matplotlib.pyplot as plt
#  plt.ioff()
import os

# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman'],
                                'monospace': ['Computer Modern Typewriter']})
rc('text', usetex=True)

# for legend
from matplotlib.font_manager import FontProperties # for legend
fontP=FontProperties()
#  fontP.set_size('large')
plt.rcParams['legend.title_fontsize'] = 'large'
fontP.set_family('serif')

print(plt.rcParams.keys())
plt.rcParams['xtick.labelsize'] = 'large'
plt.rcParams['ytick.labelsize'] = 'large'


Mpc = 3.086e24 # Mpc/cm






zlist = [ 0.1, 0.25, 0.45, 0.7, 0.9, 1.25, 1.65, 2.25 ] # list of redshift bin centers
snapshotlist = [ 61,  53, 44, 36, 31, 24, 19, 14]   # snapshot numbers corresponding to chosen redshifts
simdirs = ['galaxy-512-new-cosmo-69MPC', 'galaxy-512-new-cosmo-100MPC'] # root directories of simulations
simname = {} # dict for axis labelling
simname[simdirs[0]] = 'G69'
simname[simdirs[1]] = 'G100'



# directory where observations are
obsdir = '/home/mivkov/UZH/Masterarbeit/masterarbeit/files/observational_data/behroozi-2013-data-compilation/smf_ms/'

# list of files of observations
obsfilelist = [
    [ 'moustakas_z0.105.smf-logarithmic'   ],
    [ 'moustakas_z0.25.smf-logarithmic'    ],
    [ 'moustakas_z0.45.smf-logarithmic'    ],
    [ 'perez_gonzalez_z0.70.smf', 'moustakas_z0.725.smf-logarithmic'],
    [ 'perez_gonzalez_z0.90.smf', 'moustakas_z0.9.smf-logarithmic'], 
    [ 'mortlock_z1.0.smf'], 
    [ 'marchesini_z1.7.smf'],
    [ 'mortlock_z2.0.smf'],
 ]

obslabellist = [
    [ 'MOU '], ['MOU'], ['MOU'], ['PG', 'MOU'], ['PG', 'MOU'], ['MOR'], ['MAR'], ['MOR']
]






#==================================
def get_my_smf(i, simdir):
#==================================
    """
    Read in my SMF from smf.txt
    """

    snapdir = os.path.join(simdir, 'output_'+str(snapshotlist[i]).zfill(5))
    smffile = os.path.join(snapdir, 'smf.txt')
    infofname = os.path.join(snapdir, 'info_'+str(snapshotlist[i]).zfill(5)+'.txt')

    #------------------------------
    # first read infofile
    #------------------------------

    infofile = open(infofname)
    for j in range(9):
        infofile.readline() # skip first 9 lines

    # get expansion factor
    aline = infofile.readline()
    astring, equal, aval = aline.partition("=")
    afloat = float(aval)
    a_exp = afloat

    for j in range(5):
        infofile.readline() # skip first 9 lines

    lline = infofile.readline()
    lstring, equal, lval=lline.partition('=')
    lfloat = float(lval)
    unit_l = lfloat/Mpc

    infofile.close()

    redshift = 1.0/a_exp - 1


    #-----------------------------
    # Now read SMF
    #-----------------------------

    masses, smf = np.loadtxt(smffile, dtype='float', skiprows=1, usecols=([0, 1]), unpack=True)


    #-----------------------------
    # Process data
    #-----------------------------

    masses_plot = 0.5*(masses[1:]+masses[:-1])
    dex_divide = masses[1:]-masses[:-1]

    # Note: smf are conventionally defined over comoving volume
    vol = (unit_l*(1+redshift))**3

    smf = smf[1:] # skip first empty bin
    mask = smf > 0
    dd = dex_divide[mask]
    phi = np.log10(smf[mask]/vol/dd)
    masses = masses_plot[mask]
    
    return masses, phi




#==============
def main():
#==============

    nrows = 2
    ncols = 4

    fig = plt.figure(figsize=(10, 5))
    for i in range(nrows*ncols):
        fig.add_subplot(nrows, ncols, i+1)


    for i, ax in enumerate(fig.axes):

        # plot my simulations
        for simdir in simdirs:
            mass, phi = get_my_smf(i, simdir)
            ax.plot(mass, phi, label=simname[simdir])

        # plot observations
        for j, obs in enumerate(obsfilelist[i]):
            obsfile = os.path.join(obsdir, obs)
            print("reading obsfile", obsfile)
            mass, phi, yerrdown, yerrup = np.loadtxt(obsfile, comments='#', unpack=True, dtype=np.float) 

            ax.errorbar(mass, phi, yerr = (yerrup, yerrdown), label=obslabellist[i][j], zorder=0, capsize=2)

        axtitle = r'$z \sim {0:4.2f}$'.format(zlist[i])

        ax.legend(title=axtitle, loc='lower left', prop=fontP)
        ax.set_ylim(-6.2,-0.2)
        ax.set_xlim(5,12.5)
        ax.tick_params(axis='both', direction='in')
        ax.tick_params(axis='both', bottom=True, top=True, left=True, right=True)

        if i != 0 and i != ncols:
            ax.set_yticklabels([])

    plt.subplots_adjust(left=0.07, right=0.99, top=0.99, bottom=0.1, wspace=0., hspace=0.)
    plt.figtext(0.02, 0.5, r"$\log_{10} \Phi(M_*)$ $[($Mpc$(1+z))^{-3}$ dex$^{-1}$]", ha='center', va='center',rotation=90, size=14)
    plt.figtext(0.5, 0.015, r"$\log_{10}$ $ M_* / M_{\odot}$ ", ha='center', va='center', size=14)

    #  plt.show()
    plt.savefig('smf_both_sims-reduced.pdf', format='pdf')

    #  plt.subplots_adjust(left=0.06, right=0.99, top=0.99, bottom=0.08,wspace=0., hspace=0.)
    #  plt.figtext(0.02, 0.5, r"$\log_{10} \Phi(M_*)$ $[($Mpc$(1+z))^{-3}$ dex$^{-1}$]", ha='center', va='center',rotation=90, size=16)
    #  plt.figtext(0.5, 0.02, r"$\log_{10}$ $ M_* / M_{\odot}$ ", ha='center', va='center', size=16)
    #
    #  #  plt.show()
    #  plt.savefig('smf_both_sims-new-cosmo.pdf', format='pdf')

    return









if __name__ == "__main__":
    main()
