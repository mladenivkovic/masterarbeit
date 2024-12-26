#!/bin/bash

# The part of the makemovie script that implements the ffmpeg

#ffmpeg -f image2 -framerate 20 -pattern_type glob -i '*.png' -r 25 ../movie.mp4

#Take filename as cmd line arg.
#If none is given, take movie.mp4

if [ "$#" -gt 1 ] || [ "$#" == 0 ] ; then
    echo "movie name is: ../movie.mp4"
    filename="../movie.mp4"
else
    filename="$1"
    echo "filename is:" $filename
fi



#/users/biernack/bin/ffmpeg -r 20 -pattern_type glob -i '*.png' -c:v h264 -pix_fmt yuv420p -vf scale=720:-2  -qp 0 "$filename"
#^for dora

ffmpeg -framerate 8 -r 8 -pattern_type glob -i '*.png' -profile:v main -pix_fmt yuv420p -profile:v baseline -level 3.0 -an -c:v libx264 "$filename"

#OLD AND USELESS

# ffmpeg -r 20 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p -vf scale=800:800  -qp 0 "$filename"
#ffmpeg -framerate 20 -pattern_type glob -i '*.png' -vf "fps=30, scale=800:800" -pix_fmt yuv420p -shortest "$filename" 
#ffmpeg -framerate 20 -pattern_type glob -i '*.png' -fmt yuv420p -shortest "$filename" 

#ffmpeg -r 20 -f image2 -pattern_type glob -i '*.png' -f mp4 -q:v 0 -vcodec mpeg4 -r 20 "$filename"
