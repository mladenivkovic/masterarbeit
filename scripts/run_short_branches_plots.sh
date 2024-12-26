#!/bin/bash

# Run the plots for haloes with short main branches
# Here: For run-ntrace-100 and haloes with > 1000 particles and MBL < 10

# for halo in 5699 ; do # for fieldtest
for halo in 90361 166730 412835 477317 625108 732080 736143 753767 753787 ; do
    extract_clumpparticles.py 63 $halo
    clumppartfile=clumpparticles_00063-clump-$halo.pkl
    plot_clumpparticles.py $clumppartfile
    mkdir -p tracing_output_00063-clump-$halo
    mv clump_*"$halo"* tracing_output_00063-clump-$halo
done
