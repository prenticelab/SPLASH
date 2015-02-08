#---------- Variables ----------#
CXX = g++
CXXFLAGS = -g

MAKEFILE = Makefile

OBJS = DATA.o EVAP.o main.o

OFILES = DATA.o EVAP.o main.o
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
main.o:		DATA.h EVAP.h
