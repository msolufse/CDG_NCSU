%--------------------------------------------------------------------------
%Solve the ODE using the given parameters (pars) and initial conditions
%(Init) at the given points (xdata)
%Set ydata to global
%In this case, x and xdata are the same so don't need to interpolate ydata
%--------------------------------------------------------------------------

function [sol,rout] = model_sol(pars,x,Init)

global ODE_TOL 

options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);
sol   = ode45(@modelBasic,[x(1) x(end)],Init,options,pars);
% Evaluates solution at time-points where data are measured
sol   = deval(sol,x);

I = sol(2,:);
rout = I;

% Sensitivity wrt residual or scaled residual
%rout = (sol(2,:)' - ydata);
%rout = (sol(2,:)' - ydata)/mean(ydata);
