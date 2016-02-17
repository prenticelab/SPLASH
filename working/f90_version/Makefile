# Copyright (C) 2015, see LICENSE, Benjamin D. Stocker
# contact: benjamin.stocker@gmail.com

# Debug compilation mode with gfortran
# FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ffree-form

## Use PGF95 compiler
# FCOM=pgf95 -r8 -Mextend -g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform

# The compiler
FC = gfortran
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
# flags forall (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I/usr/include


# name of executable
PROGRAMS = test_spinup_splash test_splash

all: $(PROGRAMS)

spinup_splash.o: splash.mod.o
test_spinup_splash: splash.mod.o
test_splash.o: splash.mod.o
test_splash: splash.mod.o

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(PROGRAMS) *.o *.mod *.MOD
