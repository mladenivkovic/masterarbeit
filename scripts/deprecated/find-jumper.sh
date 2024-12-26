#!/bin/bash


#---------------------------------------------------------------
# This script finds progenitors that have found descendants
# over non-adjacent snapshots and prints those descendants.
# Only cmdline arg should be the output number of the directory
# you want to search in. This script works nicely in for loops.
#---------------------------------------------------------------


if [[ $# -ne 1 ]]; then
    echo need output number
    exit
fi

outputnr=$1


output=output_`printf "%05d" $outputnr`
# echo $output

for i in $output/mergertree*; do
    # find jumpers
    awk -v outnr=$(($outputnr-1)) \
    'NR > 1 && $3 != outnr {print "found jumper: descendant is ", $1}' $i
    # find new clumps
    # awk 'NR > 1 && $2 == 0 {print "found new clump: ", $1}' $i
    # find mergers
    # awk 'NR > 1 && $1 < 0 {print "found merger: ", $1}' $i

done

