# for debugging

# Use gfortran compiler
FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8
# debug compilation mode with gfortran
# FCOM=gfortran -ffixed-line-length-0 -fdefault-real-8 -g -fbounds-check -Wall -fbacktrace -finit-real=nan

# Use PGF95 compiler
#FCOM=pgf95 -Mextend -g -O0 -Mextend -Mbounds -Minfo -Minform=inform -Kieee -Ktrap=fp 

export FCOM

# name of executable
EXE = runstash

standard:
	$(FCOM) stash.F

.PHONY: clean
	clean:
	rm $(EXE) *.o *.mod

