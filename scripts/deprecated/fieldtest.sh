#!/bin/bash

#----------------------------------------------------------
# This script restarts a ramses simulation
# at every output. It assumes that the simulation
# terminates after the mergertrees are finished.
# To be used when you only want to check mergertree stuff
# and save time by not advancing the simulation.
#----------------------------------------------------------




DATE=`date +%F_%Hh%M`

# # rewrite back to nrestart=1
# sed -i "s/nrestart=[0-9]*/nrestart=0/" merger.nml



namelist=merger.nml
ex=ramses3d-merger


if [[ $# == 1 ]]; then

    case $1 in

        -dice)
            namelist=merger-dice.nml
            ex=ramses3d-merger-dice
        ;;

        *)
            echo "use -dice for dice runs"
            exit
    esac

fi


for out in output_*; do
    nr=$(echo ${out#output_} | sed 's/^0*//')
    echo "============================"
    echo "WORKING FOR NR", $nr
    echo "============================"
    echo; echo; echo; echo;
    sed -i "s/nrestart=[0-9]*/nrestart=${nr}/" "$namelist"
    mpiexec -n 4 "$MA"/exec/"$ex" "$namelist" 2>>"$DATE"_error.log | tee --append run"$DATE".log
    if [[ ${PIPESTATUS} -ne 0 ]]; then exit; fi;

done 

# rewrite back to nrestart=1
sed -i "s/nrestart=${nr}/nrestart=1/" "$namelist"
