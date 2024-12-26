#!/bin/bash

#---------------------------------------------------------------
# This script automatises the file transfer from daint.
# Get everything from the storage directory on daint into a
# directory on this machine.
#---------------------------------------------------------------

#----------------
# set parameters
#----------------

# root dir of where the data will be written
# rood dir = dir that contains output_* dirs
destdir=$HOME/simulation_archive/ramses/master_thesis/galaxy-512-new-cosmo-69MPC/output_00067
# destdir=$HOME/simulation_archive/ramses/master_thesis/galaxy-512-new-cosmo-100MPC/output_00067

# root dir of where to get data from
srcdir=/store/uzh/uzh5/mivkov/galaxy-512-new-cosmo/galaxy-512-69MPC/output_00067
# srcdir=/store/uzh/uzh5/mivkov/galaxy-512-new-cosmo/galaxy-512-100MPC/output_00067


[[ -d $destdir ]] || mkdir -p $destdir


for fname in correlation Pk; do
    for add in _all _sub; do
        for thresh in -0.00E+00 -1.00E+09 -1.00E+10 -1.00E+11; do
            fullf="$fname""$add""$thresh".txt

            globus-url-copy -cd -v -r -sync \
                sshftp://mivkov@gridftp.cscs.ch/$srcdir/"$fullf" \
                file://$destdir/
        done
    done
done
