# makefile for spline interpolation comparisons and timing (Lai & Kaplan, 2020)

FC=gfortran
FCFLAGS  = -Wall -pedantic

#---------------------------------------------

OBJS = parametersmod.o			\
			 utilitiesmod.o    		\
			 newsplinemod.o				\
       output_daily_final.o

.SUFFIXES: .o .f90 .f .mod

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $<

all::	newspline_test

newspline_test: $(OBJS)
	$(FC) $(FCFLAGS) -o newspline_test $(OBJS)

clean::
	-rm newspline_test *.o *.mod
