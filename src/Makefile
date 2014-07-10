F90=gfortran
#FFLAGS=-O3 # Use 3rd level of optimisation
#FFLAGS=-O2 -Wall -fdefault-real-8           # DP/ 3rd level of optimisation pedantic
FFLAGS=-O2 -fdefault-real-8 -fopenmp    # DP/ 3rd level of optimisation
   #
   #  -pg   =  profiler (gprof)
   #  -g    = debugger  (gdb)
   #

## Generic rules:
.SUFFIXES:          # get rid of that annoying Modula rule
.SUFFIXES: .f .f90 .c .o .obj .x .exe .inc .h  

#
# Rules (how to make a Blergh from a Blah
#
.f90.o:
	    $(F90) $(FFLAGS) -o $*.o -c $*.f90

#
# Special targets
#
default: run.x
#
# Targets
#
run.x: run.f90 cdata.o derivatives.o reconnection.o initial_cond.o initial.o output.o boundary.o general.o timestep.o line.o diagnostic.o normal_fluid.o tree.o Makefile 
	$(F90) $(FFLAGS) -o ./run.x run.f90 cdata.o derivatives.o reconnection.o initial_cond.o initial.o output.o boundary.o general.o timestep.o line.o diagnostic.o tree.o normal_fluid.o
#
#
clean: 
	rm -f *.o *.x *.exe *.obj *.mod
cleann: clean
	rm -f *~
pristine: cleann
	rm -f *.dat *.log *.png
#
# Dependencies
#
cdata.o: cdata.f90 ../run.in
# 
output.o: output.f90 cdata.o diagnostic.o
#
boundary.o: boundary.f90 cdata.o
#
general.o: general.f90 derivatives.o cdata.o
#
derivatives.o: derivatives.f90 cdata.o
#
initial.o: initial.f90 initial_cond.o cdata.o boundary.o output.o reconnection.o normal_fluid.o
#
initial_cond.o: initial_cond.f90 cdata.o general.o boundary.o
#
normal_fluid.o: normal_fluid.f90 cdata.o general.o 
#
tree.o: tree.f90 general.o cdata.o
#
timestep.o: timestep.f90  general.o cdata.o normal_fluid.o tree.o
#
line.o: line.f90 reconnection.o boundary.o general.o cdata.o 
#
reconnection.o: reconnection.f90 general.o cdata.o
#
diagnostic.o: diagnostic.f90  general.o cdata.o normal_fluid.o
#
