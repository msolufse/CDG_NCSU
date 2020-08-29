%--------------------------------------------------------------------------
%Solve the ODE using the given parameters (pars) and initial conditions
%(Init) at the given points (ydata)
%Set ODE tolerance and data as global
%td and xdata
%--------------------------------------------------------------------------

function [sol,rout] = model_sol(pars,data)

global ODE_TOL 

ydata = data.ydata;
xdata = data.xdata;
Init  = data.Init;

options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);
sol   = ode45(@modelBasic,xdata,Init,options,pars);
sol   = deval(sol,xdata);

I = sol(2,:)';

%rout = (I - ydata);
rout = (I - ydata)/mean(ydata);
