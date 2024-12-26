#!/bin/bash

# create plots of all orphans that are in a halo.
# usage: give me number of directory to work with

if [ $# -ne 1 ]; then
    echo I need the output directory number as cmdlinearg
    exit
fi

cwd=`pwd`
dirnr=$1
dirnr_formatted=`printf '%05d' $dirnr`

echo WARNING: MAKE SURE THE HARDCODED MASS THRESHOLD IS SAME EVERYWHERE

cd $MA/scripts

make get_clump_ID_of_orphans || exit
make get_orphan_image_data || exit

cd "$cwd"

# generate orphan_clump_IDs-XXXXX.txt files
if [ ! -f orphan_clump_IDs-$dirnr_formatted.txt ]; then
    get_clump_ID_of_orphans $dirnr | tee get_clump_ID_of_orphans.log || exit
fi

# generate orphan_main_halos-00067.txt files
if [ ! -f unique_orphan_main_halos-$dirnr_formatted.txt ]; then
    get_orphan_halos.py $dirnr || exit
fi


# generate plot for each halo, go line by line
while read line; do
    get_orphan_image_data $dirnr $line
    plot_orphans.py $dirnr $line
done < unique_orphan_main_halos-$dirnr_formatted.txt

