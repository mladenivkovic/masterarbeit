#!/bin/bash


#---------------------------------
# run ramses/utils/f90/partcenter
# for 10 most massive haloes
#---------------------------------

scr=$MA/scripts
srcdir='output_00067'

for halo in `$scr/find_mostmassive.py -n 100 -d $srcdir`; do
    x=`grep '^\s*'"$halo"'\s' "$srcdir"/halo_* | awk '{print $4}'`
    y=`grep '^\s*'"$halo"'\s' "$srcdir"/halo_* | awk '{print $5}'`
    z=`grep '^\s*'"$halo"'\s' "$srcdir"/halo_* | awk '{print $6}'`

    ~/local/ramses/utils/f90/partcenter -inp $srcdir \
                                        -out $srcdir/partcenter-output-halo-$halo.dat \
                                        -xce $x \
                                        -yce $y \
                                        -zce $z \
                                        -rma 0.02
done

# # for 5298921 - most massive halo in galaxy-512-new-cosmo-100MPC
# /home/mivkov/local/ramses/utils/f90/partcenter \
#         -inp output_00067 \
#         -out output_00067/partcenter-output-halo-5298921.dat \
#         -xce 0.695 -yce 0.015 -zce 0.74 -rma 0.02
