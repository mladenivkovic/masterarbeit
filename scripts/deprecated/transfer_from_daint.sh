#!/bin/bash

#---------------------------------------------------------------
# This script automatises the file transfer from daint.
# It only transfers necessary files for the mergertree:
# mergertree.txt*, unb_form_out*, info_*, and the logs.
# Specify in the parameters from which directory from /scratch
# you want to copy and which outputs to take.
#---------------------------------------------------------------

#----------------
# set parameters
#----------------

# destdir=/media/mivkov/DATA
destdir=$MA/results/daint

begin=41
end=41
# begin=67
# end=67
ncpu=144 # only needed for scp, you can ignore it for globus-url



cd $destdir


#=======================
# GRIDFTP
#=======================

# for i in 1 10 20 50 100 200 500 1000; do
#
#     srcdir=run-ntrace-$i

# for srcdir in run-inclusive-nosaddle run-inclusive-saddle run-exclusive-nosaddle run-exclusive-saddle; do
#

# for srcdir in galaxy-256-old-redshifts; do
# for srcdir in galaxy-512; do
# for srcdir in galaxy-512-new-redshifts; do
# for srcdir in galaxy-512-100MPC; do
for srcdir in galaxy-512-100MPC galaxy-512-69MPC; do

    [[ -d $destdir/$srcdir ]] || mkdir -p $destdir/$srcdir


    #-------------------------
    # Copy simulation files
    #-------------------------

    for out in $(seq -f "%05g" $begin $end ); do
        # for file_prefix in density_image smf part2map; do
        # for file_prefix in Pk_ correlation_; do
        # for file_prefix in galaxies_ halo_ info_ ; do
        # for file_prefix in galaxies_ halo_ clump_ info_ mergertree_ ; do
        # for file_prefix in galaxies_ ; do
        # for file_prefix in part_ unbinding_; do
        # for file_prefix in halo_ clump_ mergertree info_ part_ unbinding_; do
        # for file_prefix in density_image smf part2map halo_ galaxies_ info_  Pk_ correlation_; do
        for file_prefix in positions_ radial_galaxies_ density_image_*; do
            globus-url-copy -cd -v \
                sshftp://mivkov@nordend01.cscs.ch/scratch/snx3000/mivkov/$srcdir/output_$out/$file_prefix* \
                file://$destdir/$srcdir/output_$out/
                # sshftp://mivkov@gridftp.cscs.ch/scratch/snx3000/mivkov/$srcdir/output_$out/$file_prefix* \
                # file:///$destdir/$srcdir/output_$out/
        done
    done



    #----------------------
    # Copy logfiles
    #----------------------

    # for logfile_prefix in run slurm merger job; do
    for logfile_prefix in log_halostats; do
        globus-url-copy -cd -v\
            sshftp://mivkov@nordend01.cscs.ch/scratch/snx3000/mivkov/$srcdir/$logfile_prefix* \
            file:///$destdir/$srcdir/
    done



    #----------------------
    # Evaluate
    #----------------------
    # cd $destdir/$srcdir
    # eval_trees.py
    # # ln -n $destdir/$srcdir/eval_trees.txt $MA/results/ntracers/eval_trees-ntrace-$i.txt
    # cd $destdir


    #----------------------
    # Plot
    #----------------------

    # for halo in 677880 655230 336245 50660 291511; do
    #     treeplot.py -v -pdf -pp -lc $halo | tee log-"$halo".txt;
    # done

done







#=======================
# SCP
#=======================

# for srcdir in galaxy-256-nocalc; do
#
#     [[ -d $destdir/$srcdir ]] || mkdir -p $destdir/$srcdir
#
#
#     #-------------------------
#     # Copy simulation files
#     #-------------------------
#
#     for out in $(seq -f "%05g" "$begin" "$end" ); do
#         mkdir -p $srcdir/output_"$out"
#         # for file_prefix in halo_ clump_ mergertree info_ part_ unbinding_; do
#         for file_prefix in mergertree halo galaxies; do
#             for cpu in $(seq -f "%05g" "$ncpu" ); do
#                 fn="$srcdir"/output_"$out"/"$file_prefix"_"$out".txt"$cpu"
#                 scp daint:/scratch/snx3000/mivkov/"$fn" $destdir/$fn
#             done
#         done
#         fn="$srcdir"/output_"$out"/info_"$out".txt
#         scp daint:/scratch/snx3000/mivkov/"$fn" $destdir/"$fn"
#     done
#
#
#
#     #----------------------
#     # Copy logfiles
#     #----------------------
#
#     # for logfile_prefix in run slurm merger job; do
#     #     scp "daint:/scratch/snx3000/mivkov/$srcdir/$logfile_prefix\*" $destdir/$srcdir/
#     # done
#
#
#
#     #----------------------
#     # Evaluate
#     #----------------------
#     # cd $destdir/$srcdir
#     # eval_trees.py
#     # # ln -n $destdir/$srcdir/eval_trees.txt $MA/results/ntracers/eval_trees-ntrace-$i.txt
#     # cd $destdir
#
#
#     #----------------------
#     # Plot
#     #----------------------
#
#     # for halo in 677880 655230 336245 50660 291511; do
#     #     treeplot.py -v -pdf -pp -lc $halo | tee log-"$halo".txt;
#     # done
#
# done
#


