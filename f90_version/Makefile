## Use gfortran compiler
# FCOM=gfortran -ffixed-line-length-0
## debug compilation mode with gfortran
# FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan -ffree-form

## Use PGF95 compiler
FCOM=pgf95 -r8 -Mextend -g -O0 -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp -Mfreeform

export FCOM

# name of executable
EXE = runsplash

standard:
	$(FCOM) splash.f90 -o $(EXE)

.PHONY: clean
	clean:
	rm $(EXE) *.o *.mod

