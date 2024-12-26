#!/bin/bash


#==========================================================================
# This script compiles and calls the fortran get_smf.f03
# program correctly.
# assumes you are in the parent directory of output_XXXXX files.
# usage:
# get_halostats.sh <output_nr> <halo-Id>
#==========================================================================



#------------------------
# compile f03 program
#------------------------
MASCR=~/UZH/Masterarbeit/masterarbeit/scripts

workdir=$PWD
cd $MASCR

#local
# gfortran get_halostats.f03 -o get_halostats.o -O3 -g -fbacktrace -fbounds-check -Wall -fopenmp  -frecord-marker=4
# gfortran get_halostats.f03 -o get_halostats.o -O0 -g -fbacktrace -fbounds-check -Wall -frecord-marker=4 -fopenmp
# hydra
# gfortran get_smf.f03 -o get_smf.o -g -fbacktrace -fbounds-check -fopenmp -I/usr/include -lfftw3 -lfftw3_omp -Wall

make get_halostats

if [ $? -ne 0 ]; then
    echo COMPILATION FAILED. EXITING.
    cd $workdir
    exit
fi

cd $workdir



#--------------------------
# call fortran program
#--------------------------

get_halostats $1 $2

if [ $? -ne 0 ]; then
    echo GET_HALOSTATS.F03 FAILED. EXITING
    exit
fi



#-------------------------
# Plot halos
#-------------------------

python2 $MASCR/plot_fortran_halo.py $1 $2
if [ $? -ne 0 ]; then
    echo PLOT_FORTRAN_SMF.PY FAILED. EXITING
    exit
fi
# eog halo-density-image_output_`printf "%05d" $1`-halo-$2.png
