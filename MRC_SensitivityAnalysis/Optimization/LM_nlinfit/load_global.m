function [pars,Init,low,hi] = load_global
global ODE_TOL DIFF_INC

ODE_TOL  = 1e-8;
DIFF_INC = sqrt(ODE_TOL);

S0 = 900; I0 = 100; R0 = 0;
Init = [S0 I0];
    
gamma = 0.2;	k = 0.1;	r = 0.6;	delta = 0.15;
pars = [gamma k r delta]';


