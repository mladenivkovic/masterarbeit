#!/bin/bash

#-------------------------------------------------
# Counts how many particles of a clump are in
# the unb_form_out* files.
# Clump must be only argument on cmdline.
# Execute from within output_XXXXX directory.
#-------------------------------------------------



clumpid=$1


particles=0
for i in unb_form_out_particleoutput*; do
    echo checking file $i
    count=`grep "\s$clumpid\s" $i | wc -l`
    echo found $count particles
    particles=$(($particles + $count))
done

echo
echo Found $particles particles of clump $clumpid in total

