#!/bin/bash


#==========================================================================
# This script compiles and calls the fortran get_smf.f03
# program correctly.
# assumes you are in the parent directory of output_XXXXX files.
# usage:
# get_smf.sh
# no cmd line args!
#==========================================================================



#------------------------
# compile f03 program
#------------------------

workdir=$PWD
cd $MA/scripts

#local
gfortran get_smf_types.f03 -o get_smf_types.o -O3 -g -fbacktrace -fbounds-check -fopenmp -Wall
# hydra
# gfortran get_smf.f03 -o get_smf.o -g -fbacktrace -fbounds-check -fopenmp -I/usr/include -lfftw3 -lfftw3_omp -Wall

if [ $? -ne 0 ]; then
    echo COMPILATION FAILED. EXITING.
    cd $workdir
    exit
fi

cd $workdir



#--------------------------
# call fortran program
#--------------------------

get_smf_types.o

if [ $? -ne 0 ]; then
    echo GET_SMF.F03 FAILED. EXITING
    exit
fi



#-------------------------
# Plot SMF
#-------------------------

plot_fortran_smf.py
if [ $? -ne 0 ]; then
    echo PLOT_FORTRAN_SMF.PY FAILED. EXITING
    exit
fi
eog smf.png

