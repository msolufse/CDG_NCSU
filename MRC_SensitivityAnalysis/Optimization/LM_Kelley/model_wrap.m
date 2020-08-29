%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function [sol,rout,J] = model_wrap(pars,data)
global ALLPARS INDMAP 

tpars = ALLPARS;
tpars(INDMAP') = pars;

[sol,rout,J] = model_sol(tpars,data);
