FC=/usr/local/bin/gfortran

all: treebranch

treebranch.o: treebranch.f90 f90_tools.mod tree_pres.mod
	$(FC) -c treebranch.f90

f90_tools.o: f90_tools.f90
	$(FC) -c f90_tools.f90

f90_tools.mod: f90_tools.f90
	$(FC) -c f90_tools.f90

tree_pres.o: tree_pres.f90 f90_tools.mod
	$(FC) -c tree_pres.f90

tree_pres.mod: tree_pres.f90 f90_tools.mod
	$(FC) -c tree_pres.f90

treebranch: treebranch.o f90_tools.o tree_pres.o
	$(FC) -o treebranch treebranch.o f90_tools.o tree_pres.o

clean: 
	rm -f *.o *.mod Palpha* Pbeta* *.2d treebranch
