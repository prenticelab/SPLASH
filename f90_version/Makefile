## Use gfortran compiler
# FCOM=gfortran -ffixed-line-length-0
## debug compilation mode with gfortran
# FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan

## Use PGF95 compiler
FCOM=pgf95 -r8 -Mextend -g -O0 -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp
## debug compilation model with pgi

export FCOM

# name of executable
EXE = runstash

standard:
	$(FCOM) stash.F -o $(EXE)

.PHONY: clean
	clean:
	rm $(EXE) *.o *.mod

