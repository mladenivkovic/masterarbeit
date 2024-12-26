# This Makefile is used to compile the modules.f90 code correctly.
#
# SYNTAX:
#target: depencencies
#[tab] 	system command




#==================================
# COMPILATION TIME PARAMETERS
#==================================


EXEC=make_mergertree


#==================================
# COMPILER OPTIONS AND DIRECTIVES
#==================================


#Compiler to use
F90=gfortran -fbacktrace -g -O 
# flags for gfortran compiler
F90FLAGS= -fbounds-check -Wuninitialized -Wall -fcheck=all -ffree-line-length-none -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan




#==================================
# FILE LISTS
#==================================

#include paths. Will be followed in that order.
VPATH = ../src


TREEOBJECTS = make_mergertree.o ramses_copy.o

MODOBJECTS =  amr_parameters.o mergertree_commons.o io.o sort.o


%.o: %.f90 
	$(F90) $(F90FLAGS) -c $^ -o $@
	

$(EXEC): $(MODOBJECTS) $(TREEOBJECTS)
	$(F90) $(MODOBJECTS) $(TREEOBJECTS) -o $@ 

#==================================
#==================================

.PHONY: clean # keeps makefile from doing something with a file named clean

clean:
	rm -f *.o *.mod $(EXEC) 

