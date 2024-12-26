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
destdir=$HOME/simulation_archive/ramses/master_thesis/galaxy-512-new-cosmo-69MPC
# destdir=$HOME/simulation_archive/ramses/master_thesis/galaxy-512-new-cosmo-100MPC

# root dir of where to get data from
srcdir=/store/uzh/uzh5/mivkov/galaxy-512-new-cosmo/galaxy-512-69MPC
# srcdir=/store/uzh/uzh5/mivkov/galaxy-512-new-cosmo/galaxy-512-100MPC


[[ -d $destdir ]] || mkdir -p $destdir



globus-url-copy -cd -v -r -sync \
    sshftp://mivkov@nordend01.cscs.ch/$srcdir/* \
    file://$destdir/
