#---------- Variables ----------#
CXX = g++
CXXFLAGS = -g

MAKEFILE = Makefile

OBJS = DATA.o EVAP.o STASH.o main.o

OFILES = DATA.o EVAP.o STASH.o main.o
PFILES = main

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

#---------- Targets ----------#
all:	main

main:	$(OBJS) $(MAKEFILE)
	$(CXX) $(CXXFLAGS) $(OBJS) -o main

clean:
	rm -f $(OFILES) $(PFILES)

nice:
	rm -f $(OFILES)

#---------- Dependencies ----------#
DATA.o:		DATA.h
EVAP.o:		EVAP.h
STASH.o:	STASH.h
main.o:		DATA.h EVAP.h STASH.h
