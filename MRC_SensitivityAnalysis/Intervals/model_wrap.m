function sol = model_wrap(pars,x)
global ALLPARS INDMAP Init

tpars = ALLPARS;
tpars(INDMAP') = pars;

sol = model_sol(tpars,x,Init)';

sol = sol(:,2);