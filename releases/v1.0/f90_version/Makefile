# Copyright (C) 2015, see LICENSE
# Use PGI compiler ( pgf95 v.16.1-0)
#   The Portland Group - PGI Compilers and Groups
#   Copyright (c) 2016, NVIDIA Corporation

#---------- Variables ----------#
FC=pgf95
FCFLAGS = -r8 -Mextend -g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform
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
