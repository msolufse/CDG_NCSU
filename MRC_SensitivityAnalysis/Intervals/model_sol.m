function [sol] = model_sol(pars,x,Init)

global ODE_TOL 

options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);
sol   = ode45(@modelBasic, x, Init,options,pars);
sol   = deval(sol,x);


