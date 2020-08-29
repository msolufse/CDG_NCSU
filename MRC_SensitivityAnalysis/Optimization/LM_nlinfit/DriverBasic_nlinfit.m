function [x,histout,jachist,rout,sc] = DriverBasic_nlinfit

global ALLPARS INDMAP xdata ydata data


% Load data
load SIR_data.mat
ydata = obs;
xdata = (0:0.1:6)';

% Parameters estimated
INDMAP = [2 3 4]; 

% Get nominal parameter values, upper and lower bound for optimization
[pars,Init] = load_global;

% Create structure with data and initial conditions
data.ydata = ydata;
data.xdata = xdata;
data.Init  = Init;

%--------  Solution w/ Nominal parameter values -----------------------
% Integrate the system using ode45.m
sol = model_sol(pars,data);
Isol = sol(2,:)';

figure(1); hold on
h = plot(xdata,ydata,'*',xdata,Isol);
set(h(1),'linewidth',4);
set(h(2),'linewidth',8);
set(h,'Markersize',7);
set(gca,'Fontsize',24);
xlabel('Time')
ylabel('Number of infections')
grid on

%----------------------------------nlinfit---------------------------------
pars(2:end) = pars(2:end).*rand(size(pars(2:end)));
ALLPARS = pars;

w = ones(1,length(xdata))';
%use if 
%r = (ymodel-ydata)/mean(ydata) 
%w = ones(1,length(xdata))'/(mean(ydata))^2; 

[xopt,r,J,COV,MSe] = nlinfit(xdata,ydata,@model_wrap,pars(INDMAP),'Weights',w);

parsLM2 = ALLPARS;
parsLM2(INDMAP) = xopt

%Integrate the system using ode45.m
sol = model_sol(parsLM2,data);
Isol = sol(2,:)';

figure(1); hold on
hhhh = plot(xdata,Isol);
set(hhhh,'linewidth',4);
%legend([h(1) h(2) hhhh], 'Data', 'Nominal', 'Nelder-Mead','LM (nlinfit)')
print -depsc2 opt_nonlinfit.eps

save optResults.mat parsLM2 r J COV MSe