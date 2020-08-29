%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function sol = model_wrap(pars)
global ALLPARS INDMAP data

tpars = ALLPARS;
tpars(INDMAP') = pars;

sol = model_sol(tpars,data)';
