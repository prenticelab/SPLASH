# Copyright (C) 2016 Prentice Lab
# Use GNU C++ compiler (g++ v.4.8.2)
#   GNU Compiler Collection
#   Copyright (C) Free Software Foundation, Inc.

#---------- Variables ----------#
CXX = g++
CXXFLAGS = -g

MAKEFILE = Makefile

OBJS = DATA.o SOLAR.o EVAP.o SPLASH.o global.o main.o

OFILES = DATA.o SOLAR.o EVAP.o SPLASH.o global.o main.o
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
SOLAR.o:    global.h SOLAR.h
EVAP.o:		global.h SOLAR.h EVAP.h
SPLASH.o:	global.h DATA.h SOLAR.h EVAP.h SPLASH.h
global.o:   global.h
main.o:		DATA.h SOLAR.h EVAP.h SPLASH.h
