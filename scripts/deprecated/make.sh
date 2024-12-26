#!/bin/bash

FLAGS="-frecord-marker=4 -ffree-line-length-none -fbacktrace -g -O -fbounds-check -Wuninitialized -Wall -x f95-cpp-input -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan"

mpif90 $FLAGS read_progenitors.f90 -o read_progenitors.o
