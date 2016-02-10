## Use gfortran compiler
# FCOM=gfortran -ffixed-line-length-0
## debug compilation mode with gfortran
# FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ffree-form

## Use PGF95 compiler
FCOM=pgf95 -r8 -Mextend -g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform

export FCOM

# name of executable
EXE = runsplash

standard:
	$(FCOM) -c splash.mod.f90
	$(FCOM) splash.mod.o spinup_splash.f90 -o $(EXE)	

test:
	$(FCOM) -c splash.mod.f90
	$(FCOM) splash.mod.o test_splash.f90 -o $(EXE)

.PHONY: clean
clean:
	rm $(EXE) *.o *.mod

