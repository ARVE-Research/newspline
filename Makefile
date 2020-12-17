# makefile for spline interpolation comparisons and timing (Lai & Kaplan, 2020)

FC=gfortran
FCFLAGS  = -ffree-form -ffree-line-length-none

#---------------------------------------------

OBJS = parametersmod.o			\
			 utilitiesmod.o    		\
			 newsplinemod.o				\
       output_daily_final.o


.SUFFIXES: .o .f90 .f .mod

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $<

all::	newspline

newspline: $(OBJS)
	$(FC) $(FCFLAGS) -o newspline $(OBJS)

clean::
	-rm newspline *.o *.mod
