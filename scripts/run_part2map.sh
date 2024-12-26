#!/bin/bash
#------------------------------------------------------
# A script to run part2map and plot images directly.
# Make sure you compiled part2map first!
# (in ramses/utils/f90/ , use Makefile)
#
#
# usage: run_part2map.sh <optional: nc>
# you need to be in directory that contains output_XXXXX
# directories
#------------------------------------------------------

nc=2000

# set the path to the executable here
ex=$HOME/local/ramses/utils/f90/part2map

# for out in output_*; do
for out in output_00041; do 

    $ex -inp $out -out $out/part2map.dat -nx $nc -ny $nc -fil bin -per T
    python3 plot_part2map.py $out
    # eog particleplot_${out#output_}.png

done
