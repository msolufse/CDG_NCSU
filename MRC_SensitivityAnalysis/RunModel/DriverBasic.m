%--------------------------------------------------------------------------
%Computes the sensitivity matrix dy/dpars.
%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set 
%pars = log(pars) in DriverBasic_sens
%Plots ranked sensitivities 
%--------------------------------------------------------------------------

function DriverBasic

global ydata xdata ODE_TOL

load SIR_data.mat %obs

[pars,Init] = load_global;
xdata = (0:0.1:6)';
ydata = obs;

size(obs)
size(xdata)

options = odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);
sol   = ode45(@modelBasic,[xdata(1) xdata(end)],Init,options,pars);

time = sol.x;
S = sol.y(1,:);
I = sol.y(2,:);

figure(1);clf; 
h=plot(time,S);
set(h,'linewidth',4);
set(gca,'Fontsize',24);
xlabel('time (days)');
ylabel('Susceptible Individuals');
grid on;
print -depsc2 S.eps

figure(2); clf;
h=plot(time,I,xdata,ydata,'r*');
set(h,'Linewidth',4);
set(gca,'Fontsize',24);
xlabel('time (days)');
ylabel('Infected individuals');
grid on;
print -depsc2 I.eps