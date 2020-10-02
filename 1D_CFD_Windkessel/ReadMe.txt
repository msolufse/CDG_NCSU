Intro:

This folder contains a combined C++ and Fortran code for a single vessel model of pulse wave propagation. The example code uses 4-element Windkessel model that is reducible to a three element Windkessel model. A generic inflow waveform is used to drive the system. A Makefile is used to link all the modules in the right order. Current version of the Makefile is updated to work on an iMac (Make sure to modify the compiler versions and/or path links specific to your machine).



Files Description

1. Makefile: Run it by typing "make" command in the command window to create an executable “sor06”. Google “Visual Studio For Windows” for help making on the Windows. There might be many other options available.
2. main.m: Run this file in the MATLAB After creating the executable. Change geometric and hemodynamic parameter in this film to create a new example

C++ Files

3. sor06.h:  Header file containing global parameters
4. sor06.C:  Defines the network and calls the solver and the prints the data in vessel specific files.
5. arteries.h: Header file declaring the class “Tube” and all its objects and functions, to be specified in arteries.C
6. arteries.C: Contains the main computational code, which uses Lax-Wenderoff scheme and a linear pressure-area relation. It also calls tools.C, tools.h (for root finding at bifurcations) and .f90 routines to compute impedance for the outflow condition.
7. tools.C: Includes numerical algebraic tools for finding roots for the system of equations
8. tools.h: Used by tools.C and arteries.C

Fortran Files

9.  root_imp.f90: Computes input impedance using the Windkessel model at the terminal point of large arteries
10. impedance_init_sub.f90 and
11. impedance_sub.f90: Initiate and organize the impedance calculation
12. f90_tools.f90: Contains FFT and IFFT routines

Input data

13. Qin_8192.dat: Generic inflow data created from InFlow.m

Note: The number of time steps in the input data should be an integer power of 2 (To meet the requirements of FFT in f90_tools.f90)

Other MATLAB files

14. Example.m: Plots the simulated data as 2D and 3D graphs.
15. gnuplot.m: Rearrange the block matrix output
16. InFLow.m: Generates an inflow profile and saves it as Qin_8192.dat.
17. ND_Par.m: Computes and non-dimensionalize vessel stiffness and the Windkessel parameters.

Output files

Code will generate a new data files with extension *.2d



