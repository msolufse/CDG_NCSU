function [x,histout,jachist,rout,sc] = DriverBasic_fminsearch

global ALLPARS INDMAP data

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

figure(1)
hold on
h = plot(xdata,ydata,'*',xdata,Isol);
set(h(1),'linewidth',4);
set(h(2),'linewidth',8);
set(h,'Markersize',7);
set(gca,'Fontsize',24);
xlabel('Time')
ylabel('Number of infections')
grid on

%----------------------- Nelder-Mead (fminsearch) -----------------------
pars(2:end) = pars(2:end).*rand(size(pars(2:end)));
ALLPARS = pars;

opts = optimset('MaxIter',5000,'MaxFunEvals',5000);
[xopt, ~, ~, ~] = fminsearch(@model_fmin,pars(INDMAP),opts);

parsNM = ALLPARS;
parsNM(INDMAP) = xopt

% Integrate the system using ode45.m
sol = model_sol(parsNM,data);
Iopt = sol(2,:)';

figure(1); hold on
hh = plot(xdata,Iopt);
set(hh,'linewidth',3);
legend([h(1) h(2) hh], 'Data', 'Nominal', 'Nelder-Mead');
print -depsc2 FminsearchOpt.eps

save optResults.mat parsNM Iopt 
