#!/usr/bin/python3


#===================================================
# This script plots the results of eval_tree.py
# tailor-made for the presentation.
# saves them as a pdf
#===================================================

import numpy as np
import matplotlib.pyplot as plt
# use LaTeX text
from matplotlib import rc
rc('font', **{'family':'serif',
    'serif':['Computer Modern Roman'],
    'monospace': ['Computer Modern Typewriter'], 'size':16})
rc('text', usetex=True)



#=========================================================
def get_line_to_array(linestring, maxcounter, toint=False):
#=========================================================
    """
    Translate string line separated by commas to numpy array of integers (if toint=True) or floats (if toint=False)
    simultaneously find highest value in array
    """

    import numpy as np

    linestring = linestring.strip()
    outlist = linestring.split(",")

    if toint:
        for i, val in enumerate(outlist):
            outlist[i] = int(val.strip())
            maxcounter = max(maxcounter, outlist[i])
        out = np.array(outlist, dtype='int')

    else:
        for i, val in enumerate(outlist):
            outlist[i] = float(val.strip())
            maxcounter = max(maxcounter, outlist[i])

        out = np.array(outlist, dtype='float')

    return out, maxcounter





#=================================
def save_figures(figname, fig):
#=================================
    """
    Saves figure fig as figname in .png and .tex format.
    figh, figw: string for height/width of tex figure
    """


    fig.tight_layout()
    fig.savefig(figname, dpi=300, form='pdf')
    print("Saved", figname)

    return




#===================
def main():
#===================


    #================================
    # Select files
    #================================

    #--------------------
    # Debug
    #--------------------
    #  allfiles = ['eval_trees.txt']
    #  labelnames = ['label']


    #--------------------
    # Actually in use
    #--------------------
    #  ntrace = [1, 10, 50, 100, 200, 500, 1000]
    #  allfiles = ['eval_trees-ntrace-'+str(i)+'.txt' for i in ntrace]
    #  labelnames = [r'$n_{mb}='+str(i)+'$' for i in ntrace]

    allfiles = [ 'eval_trees-inc-nosaddle.txt' ] #,  'eval_trees-exc-nosaddle.txt' ]
    labelnames =[ 'inclusive no saddle' ] #, 'exclusive no saddle' ]

    #  linestyle = ['-', ':', '--', '-.', '-', ':', '--', '-.']
    linestyle = ['-']*10
    #  linestyle = [':']*8
    #  linestyle = ['--']*8
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#e377c2', '#bcbd22', '#17becf']

    # columns for legends
    if len(allfiles) <= 4:
        ncols = 1
    else:
        ncols = 2



    #==========================
    # Set up
    #==========================


    fig1 = plt.figure(1, figsize=(8,10))
    ax1 = fig1.add_subplot(211)
    ax2 = fig1.add_subplot(212, sharex=ax1)


    fig2 = plt.figure(2, figsize=(8,10))
    ax5 = fig2.add_subplot(211)
    ax6 = fig2.add_subplot(212, sharex=ax5)


    #  fig3 = plt.figure(3, figsize=(12,6))
    #  ax9 = fig3.add_subplot(121)
    #  ax10 = fig3.add_subplot(122)
    #
    #  fig4 = plt.figure(4, figsize=(6,6))
    #  ax11 = fig4.add_subplot(111)

    axgroup1 = [ax1, ax2]
    axgroup2 = [ax5, ax6]
    maxlen = 0
    maxcount1 = 0
    maxnbranch = 0
    maxcount2 = 0
    mincount1 = 1
    mincount2 = 1
#      mgmax = 0
    #  mgmin = 1
    #  mgcountmin=1
    #  mgcountmax=0
    #  mfmax = 0
    #  mfmin = 0
    #  mfcountmin = 1
    #  mfcountmax = 0
    #  displmax = 0
    #  displmin = 1
    #  dcountmax = 0
    #  dcountmin = 1
#




    #=============================
    # Read in data
    #=============================

    for f,srcfname in enumerate(allfiles):

        srcfile = open(srcfname, 'r')

        # skip 12 lines first
        for i in range(12):
            srcfile.readline()

        #----------------
        # read redshift
        #----------------
        redshift = srcfile.readline()
        redshift = redshift.strip()
        redshift = redshift.split(",")
        for i,z in enumerate(redshift):
            redshift[i] = float(z.strip())
        redshift = np.array(redshift)

        #----------------
        # read bins
        #----------------
        bins = srcfile.readline()
        bins = bins.strip()
        bins = bins.split(",")
        for i,b in enumerate(bins):
            bins[i] = int(b.strip())

        bins = np.array(bins)




        #------------------------------
        # length of main branch
        #------------------------------

        axcount = 0

        for i in range(1,bins.shape[0]):

            mblen_counts, maxcount1 = get_line_to_array(srcfile.readline(), maxcount1, False)
            mincount1 = min(mincount1, mblen_counts[mblen_counts>0].min())
            # replace zeroes with value orders of magnitude smaller
            mblen_counts[mblen_counts==0] = mincount1*1e-6

            mblens, maxlen = get_line_to_array(srcfile.readline(), maxlen, False)

            if i>1 and i < bins.shape[0]-1:
                continue

            axgroup1[axcount].semilogy(mblens[1:], mblen_counts[1:], lw=4, c=colors[f], ls=linestyle[f], label=labelnames[f], alpha=0.6, zorder=f+1)
            axcount += 1



        #----------------------
        # number of branches
        #----------------------

        axcount = 0

        for i in range(1,bins.shape[0]):

            nbranch_counts, maxcount2 = get_line_to_array(srcfile.readline(), maxcount2, False)
            mincount2 = min(mincount2, nbranch_counts[nbranch_counts>0].min())
            # replace zeroes with value orders of magnitude smaller
            nbranch_counts[nbranch_counts==0] = mincount2*1e-6

            nbranches, maxnbranch = get_line_to_array(srcfile.readline(), maxnbranch, toint=False)

            if i>1 and i < bins.shape[0]-1:
                continue

            axgroup2[axcount].loglog(nbranches, nbranch_counts, label=labelnames[f], lw=4, c=colors[f], ls=linestyle[f], alpha=0.7)
            axcount += 1





        #-------------------------
        # Mass growth
        #-------------------------

        #  mass_growth_counts, mgcountmax = get_line_to_array(srcfile.readline(), mgcountmax, False)
        #  mgcountmin = min(mass_growth_counts[mass_growth_counts>0].min(), mgcountmin)
        #  # replace zeroes with value orders of magnitude smaller
        #  mass_growth_counts[mass_growth_counts==0] = mgcountmin*1e-6
        #
        #  mass_growth, mgmax = get_line_to_array(srcfile.readline(), mgmax, False)
        #  mgmin = min(mgmin, mass_growth.min())
        #
        #  ax9.semilogy(mass_growth, mass_growth_counts, label=labelnames[f], lw=2, c=colors[f], ls=linestyle[f], alpha=0.7)
        #

        #----------------------
        # mass fluctuations
        #----------------------

        #  mass_fluct_counts, mfcountmax = get_line_to_array(srcfile.readline(), mfcountmax, False)
        #  mfcountmin = min(mass_fluct_counts[mass_fluct_counts>0].min(), mfcountmin)
        #
        #
        #  # replace zeroes with value orders of magnitude smaller
        #  mass_fluct_counts[mass_fluct_counts==0] = mfcountmin*1e-6
        #
        #
        #  mass_flucts, mfmax = get_line_to_array(srcfile.readline(), mfmax, False)
        #  mfmin = min(mass_flucts.min(), mfmin)

        #  ax10.semilogy(mass_flucts, mass_fluct_counts, label=labelnames[f], lw=2, c=colors[f], ls=linestyle[f], alpha=0.7)





        #----------------------
        # Displacements
        #----------------------

        #  displacement_counts, dcountmax = get_line_to_array(srcfile.readline(), dcountmax, False)
        #  dcountmin = min(displacement_counts[displacement_counts>0].min(), dcountmin)
        #  # replace zeroes with value orders of magnitude smaller
        #  displacement_counts[displacement_counts==0] = dcountmin * 1e-6
        #
        #
        #
        #  displacements, junk = get_line_to_array(srcfile.readline(), 0, False)
        #  displmax = max(displacements[displacement_counts>0].max(), displmax)
        #
        #  ax11.semilogy(displacements, displacement_counts, label=labelnames[f], lw=2, c=colors[f], ls=linestyle[f], alpha=0.7)



        srcfile.close()











    #==================================
    # Tweak plots and save figures
    #==================================



    # round counts
    maxlen = int(maxlen + 0.5)
    maxnbranch = int(maxnbranch + 0.5)

    #-------------------------
    # Main branch length
    #-------------------------


    xmax = maxlen+1
    xstep = 10
    ymax = maxcount1 + maxcount1//10
    ymin = mincount1-mincount1//10
    ylabelmax = int(np.log10(ymax)+0.5)+1
    ylabelmin = int(np.log10(ymin))

    for i,ax in enumerate(axgroup1):
        ax.grid()
        ax.set_xlim(0.5, maxlen+0.5)
        ax.set_xticks([i for i in range(1, xmax, xstep)])
        ax.set_xticklabels([i for i in range(1, xmax, xstep)])
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel(r'$N/N_{tot}$')
        ax.set_yticks([10**i for i in range(ylabelmin, ylabelmax)])
        ax.set_yticklabels([r"$10^{"+str(i)+"}$" for i in range(ylabelmin, ylabelmax)])
        axtwin = ax.twinx()
        axtwin.set_yticks([])
        if i == 0:
            axtwin.set_ylabel(r"$<${0:4d} particles".format(bins[i+1]-1))
        elif i == 1:
            axtwin.set_ylabel(r"$>${0:5d} particles".format(1000))
            #  ax.legend(loc='upper left', ncol=ncols)
        #  else:
        #      axtwin.set_ylabel("{0:4d} - {1:4d} particles".format(bins[i]-1, bins[i+1]-1))
        #
    ax = ax1.twiny()
    ax.set_xlim(0.5, maxlen+0.5)
    ax.set_xticks([i for i in range(1, xmax, xstep)])
    ax.set_xticklabels(["%.2f" % abs(redshift[i]) for i in range(1, xmax, xstep)])
    ax.set_xlabel(r'redshift $z$')
    ax2.set_xlabel(r'lenght of main branch')


    figname="presentation-length_of_main_branch.pdf"
    save_figures(figname, fig1)




    #-------------------------
    # number of branches
    #-------------------------

    xmax = maxnbranch+maxnbranch/10
    xlabelmax = int(np.log10(xmax)+0.5)+1
    ymax = maxcount2*1.1
    ymin = mincount2*1.1
    ylabelmax = int(np.log10(ymax)+0.5)
    ylabelmin = int(np.log10(ymin))

    for i,ax in enumerate(axgroup2):
        ax.grid()
        ax.set_xlim(1, xmax)
        ax.set_xticks([10**i for i in range(xlabelmax)])
        ax.set_xticklabels([r"$10^{"+str(i)+"}$" for i in range(xlabelmax)])
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel(r'$N/N_{tot}$')
        ax.set_yticks([10**i for i in range(ylabelmin, ylabelmax)])
        ax.set_yticklabels([r"$10^{"+str(i)+"}$" for i in range(ylabelmin, ylabelmax)])
        axtwin = ax.twinx()
        axtwin.set_yticks([])
        if i == 0:
            axtwin.set_ylabel(r"$<${0:4d} particles".format(bins[i+1]-1))
            #  ax.legend(ncol=ncols)
        elif i == 1:
            axtwin.set_ylabel(r"$>${0:5d} particles".format(1000))
            ax.set_xlabel('number of branches')
        #  else:
        #      axtwin.set_ylabel("{0:4d} - {1:4d} particles".format(bins[i]-1, bins[i+1]-1))
        #
    figname = "presentation-number_of_branches.pdf"
    save_figures(figname, fig2)






    #-------------------------
    # Mass growth and flucts
    #-------------------------

    #  ax9.set_title("Logarithmic Mass Growth")
    #  ax9.set_xlim(mgmin*1.1, mgmax*1.1)
    #  ax9.set_xticks(np.linspace(-1,1, 11))
    #  ax9.set_xticklabels(["%.2f" % i for i in np.linspace(-1,1, 11)])
    #  ax9.set_xlabel(r'$\beta_M = \frac{2}{\pi}\arctan\frac{d \log M}{d \log t}$')
    #  ylabelmin = int(np.log10(mgcountmin))
    #  ylabelmax = int(np.log10(mgcountmax)+0.5)
    #  ax9.set_ylim(mgcountmin/2, mgcountmax*2)
    #  ax9.set_yticks([10.0**i for i in range(ylabelmin, ylabelmax, 1)])
    #  ax9.set_yticklabels([r"$10^{"+str(i)+"}$" for i in range(ylabelmin, ylabelmax, 1)])
    #  ax9.set_ylabel(r'$N/N_{tot}$')
    #  ax9.legend(loc='lower right', ncol=ncols)
    #  ax9.grid()
    #
    #  ax10.set_title("Mass Growth Fluctuations")
    #  ax10.set_xlim(mfmin*1.1, mfmax*1.1)
    #  xtickstart = int(mfmin*10-0.5)/10
    #  xtickend = int(mfmax*10+0.5)/10
    #  xticks = np.linspace(xtickstart, xtickend, int((xtickend-xtickstart)/0.2+0.5)+1)
    #  ax10.set_xticks(xticks)
    #  ax10.set_xticklabels(["%.2f" % x for x in xticks])
    #  ax10.set_xlabel(r'$\xi_M = \frac{\beta_M(k, k+1) - \beta_M(k-1, k)}{2}$')
    #  ylabelmin = int(np.log10(mfcountmin))
    #  ylabelmax = int(np.log10(mfcountmax)+0.5)
    #  ax10.set_ylim(mfcountmin/2, mfcountmax*2)
    #  ax10.set_yticks([10.0**i for i in range(ylabelmin, ylabelmax, 1)])
    #  ax10.set_yticklabels([r"$10^{"+str(i)+"}$" for i in range(ylabelmin, ylabelmax, 1)])
    #  ax10.set_ylabel(r'$N/N_{tot}$')
    #  ax10.legend(loc='upper right', ncol=ncols)
    #  ax10.grid()



    # draw 90% lines
    #  ymin, ymax = ax9.get_ylim()
    #  xl, xr = get_99_percent_position(mass_growth_counts, mass_growth)
    #  ax9.plot([xl,xl], (ymin, ymax), 'r')
    #  ax9.plot([xr,xr], (ymin, ymax), 'r')
    #
    #
    #  ymin, ymax = ax10.get_ylim()
    #  xl, xr = get_99_percent_position(mass_fluct_counts, mass_flucts)
    #  ax10.plot([xl,xl], (ymin, ymax), 'r')
    #  ax10.plot([xr,xr], (ymin, ymax), 'r')

    #  plt.sca(ax9)
    #  plt.xscale('symlog')


    #  figname="mass_growth_and_fluctuation.png"
    #  figh = "12cm"
    #  figw = "18cm"
    #  save_figures(figname, fig3, figh, figw)




    #----------------------
    # Displacements
    #----------------------

    #  ax11.legend(ncol=ncols)
    #  ax11.grid()
    #  ax11.set_ylabel(r'$N/N_{tot}$')
    #  ax11.set_xlabel(r'$\Delta r$')
    #  ax11.set_xlim(-0.1, displmax+displmax/10)
    #  ax11.set_title("Displacements")
    #
    #  xtickstart = 0
    #  xtickend = int(displmax*10)/10
    #  xticks = np.linspace(xtickstart, xtickend, int(xtickend/0.3+0.5)+1)
    #  ax11.set_xticks(xticks)
    #  ax11.set_xticklabels(["%.1f" % x for x in xticks])
    #
    #  ylabelmin = int(np.log10(dcountmin))
    #  ylabelmax = int(np.log10(dcountmax)+0.5)
    #  ax11.set_ylim(dcountmin/2, dcountmax*2)
    #  ax11.set_yticks([10.0**i for i in range(ylabelmin, ylabelmax, 1)])
    #  ax11.set_yticklabels([r"$10^{"+str(i)+"}$" for i in range(ylabelmin, ylabelmax, 1)])
    #
    #  figname="displacements.png"
    #  figh = "12cm"
    #  figw = "13cm"
    #  save_figures(figname, fig4, figh, figw, special_case='displ')


if __name__ == "__main__":
    main()
