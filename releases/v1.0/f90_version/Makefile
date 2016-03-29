# Copyright (C) 2016 Prentice Lab
# Use PGI compiler ( pgf95 v.16.1-0)
#   The Portland Group - PGI Compilers and Groups
#   Copyright (c) 2016, NVIDIA Corporation

#---------- Variables ----------#
## >>> un-comment the following 2 lines for compilation with PGI compiler
# FC=pgf95
# FCFLAGS = -r8 -Mextend -g -O0 -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform
## <<<

## >>> un-comment the follwing 2 lines for compilation with (free) gfortran compiler
FC=gfortran
FCFLAGS = -ffree-line-length-0 -ffree-form -fdefault-real-8 -fdefault-integer-8
## <<<

PROGRAMS = spinup_splash test_splash

all: $(PROGRAMS)

#---------- Dependencies ----------#
spinup_splash.o: splash.mod.o
spinup_splash: splash.mod.o
test_splash.o: splash.mod.o
test_splash: splash.mod.o

#---------- Targets ----------#
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean
clean:
	rm -f $(PROGRAMS) *.o *.mod *.MOD
