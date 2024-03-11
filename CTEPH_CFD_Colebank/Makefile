#*************************************************************************#
#                                                                         #
#  Program: Makefile                                                      #
#  Version: 2.0                                                           #
#  By: Mette Olufsen                                                      #
#  Date: 14. Jan. 1997                                                    #
#                                                                         # 
#  A makefile that ensures that all modules are linked together in the    #
#  right order.                                                           #
#*************************************************************************#

# $Id: Makefile,v 1.9 2014/10/17 12:15:01 heine Exp $
# Last modified on October 23, 2014 by M. Umar Qureshi

# FOR NCSU CLUSTER
#CXX=/usr/bin/g++
# FOR WINDOWS LAPTOP
#CXX=/mingw64/bin/g++
CXX = g++

CXXFLAGS=-O2 -g -Wall -D_REENTRANT

# FOR CLUSTER
#FC=/usr/bin/gfortran

# FOR WINDOWS LAPTOP
#FC=/mingw64/bin/gfortran
FC = gfortran
FFLAGS=-O2 -g -Wall
FLIBS=-lgfortran -lquadmath

LIBS=$(FLIBS) -lm

LDFLAGS=-O2

OBJS1=tools.o junction.o sor06.o arteries.o
OBJS2=impedance_sub.o impedance_init_sub.o root_imp.o f90_tools.o

MAIN=sor06

all: $(MAIN)

$(MAIN): $(OBJS1) $(OBJS2) 
	$(CXX) -o $(MAIN) $(LDFLAGS) $(OBJS1) $(OBJS2) $(LIBS)
	
sor06.o: sor06.c sor06.h
	$(CXX) -c $(CXXFLAGS) sor06.c

junction.o: junction.c junction.h arteries.h tools.h
	$(CXX) -c $(CXXFLAGS) junction.c
	
arteries.o: arteries.c arteries.h junction.h tools.h sor06.h
	$(CXX) -c $(CXXFLAGS) arteries.c
	
tools.o: tools.c tools.h
	$(CXX) -c $(CXXFLAGS) tools.c
		
root_imp.o: root_imp.f90 f90_tools.o
	$(FC) -c $(FFLAGS) root_imp.f90
	
f90_tools.o: f90_tools.f90
	$(FC) -c $(FFLAGS) f90_tools.f90
	
impedance_sub.o: impedance_sub.f90 f90_tools.o root_imp.o
	$(FC) -c $(FFLAGS) impedance_sub.f90
	
impedance_init_sub.o: impedance_init_sub.f90 f90_tools.o root_imp.o
	$(FC) -c $(FFLAGS) impedance_init_sub.f90
	
clean:
	-rm -f *.o *.mod *.2d Zhat*
	
veryclean: clean
	-rm $(MAIN) a.out *~ STnetwork*
