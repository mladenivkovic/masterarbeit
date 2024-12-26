#!/usr/bin/python3


#=================================================
# Read in and print out Fortran unformatted file.
# Assumes file contains integers.
# Usage:
# readfort.py <file>
# options:
#   -r  read reals
#=================================================



from scipy.io import FortranFile
from sys import argv


filename = argv[-1]

f = FortranFile(filename, 'r')
data = f.read_ints()
for i in data:
    print(i)
