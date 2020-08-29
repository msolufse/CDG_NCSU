function [x,histout,jachist,rout,sc] = DriverBasic_LM_Kelley

global ALLPARS INDMAP xdata ydata Init

% Load data
load SIR_data.mat
ydata = obs;
xdata = (0:0.1:6)';

% Parameters estimated
INDMAP = [2 3 4]; 

% Get nominal parameter values, upper and lower bound for optimization
[pars,Init,low,hi] = load_global;

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

%----------------------------Levenberg-Marquardt---------------------------
pars(2:end) = pars(2:end).*rand(size(pars(2:end)));
ALLPARS = pars;

%set upper and lower bounds
optx   = pars(INDMAP);
opthi  = hi(INDMAP);
optlow = low(INDMAP);
 
% optimization information
maxiter = 40;                   % max number of iterations        
mode    = 2;                    % Performs Levenberg-Marquart optimization
nu0     = 2.d-1;                % Regularization parameter
[xopt, histout, costdata, jachist, xhist, rout, sc] = ...
     newlsq_v2(optx,'opt_wrap',1.d-4,maxiter,mode,nu0,...
               opthi,optlow,data);

parsLM = ALLPARS;                    % only optimized parameters
parsLM(INDMAP) = xopt               % all parameters, including optimized one

% Integrate the system using ode45.m
sol = model_sol(parsLM,data);
Isol = sol(2,:)';

figure(1); hold on
hhh = plot(xdata,Isol);
set(hhh,'linewidth',4);
set(hhh,'Markersize',4);
%legend([h(1) h(2) hhhh], 'Data', 'Nominal','Levenberg-Marquardt');
print -depsc2 OptSol_LM.eps

save optResultsWeighted.mat parsLM histout jachist rout sc