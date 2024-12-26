#!/bin/bash


#==========================================================================
# This scripts compiles and calls the fortran get_correlation.f03
# program correctly.
# Assumes you are in the parent directory of output_XXXXX files.
# usage:
# get_correlation.sh <output-nr> <ncells> <case>
#   output-nr:  output number of snapshot to work with
#   ncells:     optional; number of cells to divide simulation box in
#               default: 256
#==========================================================================

# gfortran get_correlation.f03 -o get_correlation.o -g -fbacktrace -fbounds-check -fopenmp -I/usr/local/include -I/$FFTW_ROOT/include -lfftw3 -lfftw3_omp -Wall

nc=256



if [ $# -lt 3 ]; then
    case $1 in
        ''|*[!0-9]*) echo argument 1 "(" $1 ")" is not an integer?; exit;;
    esac

    output=`printf "output_%05d" "$1"`
    if [ -d "$output" ]; then
        echo "working for directory" $output
    else
        echo "didnt find directory" $output
        exit
    fi

    if [ $# -gt 1 ]; then
        case $2 in
            ''|*[!0-9]*) echo argument 2 "(" $2 ")" is not an integer?; exit;;
        esac
        nc=$2
    fi

else
    echo "I need one argument: the output number to work with"
    echo "optional second argument: ncells"
    exit
fi





#------------------------
# compile f90 program
#------------------------

workdir=$PWD

# for srcfile in get_correlation.f03 get_correlation-singledir-backup.f03; do
cd $MA/scripts

# gfortran eval_galaxies.f03 -o eval_galaxies.o -g -fbacktrace -fbounds-check -I/usr/local/include -lfftw3 -Wall
# gfortran $srcfile -o get_correlation.o -g -fbacktrace -fbounds-check -fopenmp -I/usr/local/include -lfftw3 -lfftw3_omp -Wall
# gfortran get_correlation.f03 -o get_correlation.o -g -fbacktrace -fbounds-check -fopenmp -I/usr/local/include -I/$FFTW_ROOT/include -lfftw3 -lfftw3_omp -Wall

make get_correlation

if [ $? -ne 0 ]; then
    echo COMPILATION FAILED. EXITING.
    cd $workdir
    exit
fi

cd $workdir




#--------------------------
# call fortran program
#--------------------------

# get_correlation $1 $nc

if [ $? -ne 0 ]; then
    echo GET_CORRELATION.F03 FAILED. EXITING
    exit
fi



#----------------------------
# Plot Correlations
#----------------------------

plot_fortran_correlation-with-thresholds.py $output
if [ $? -ne 0 ]; then
    echo PLOT_FORTRAN_CORRELATION.PY FAILED. EXITING
    exit
fi
# eog correlation_${output#output_}.png
eog correlations-new-cosmo.png

# done

# plot_fortran_galaxies.py
# eog galaxy_density_projection_plot.png
