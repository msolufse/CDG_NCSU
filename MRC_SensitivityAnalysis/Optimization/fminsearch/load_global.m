function [pars,Init] = load_global
global ODE_TOL 

ODE_TOL  = 1e-8;

S0 = 900; I0 = 100; R0 = 0;
Init = [S0 I0];
    
gamma = 0.2;	k = 0.1;	r = 0.6;	delta = 0.15;
pars = [gamma k r delta]';


