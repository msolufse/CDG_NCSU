%=============================================================*
%                                                             
% This is a driver file that will run the systemic fluids model from c++ and
% fortran by passing parameter values needed by the model.
%
%=============================================================*
function run_1D(sample_number)

%%
%% Ensure the make file is compiled
!make clean
!make
% make the file executable 
! chmod +x sor06

%% Define the connectivity matrix

% Converging full network
% writeconn  	 = [0 1 2 0
%                 1 3 4 0
%                 2 5 6 0
%                 3 0 0 0
%                 4 0 0 0
%                 5 0 0 0
%                 6 0 0 0];
% % %             
% terminal = [3 4 5 6];
% 
% 
% dimensions = [3.58  1.27  1.27;
%               6.24  1.19  1.19;
%               5.58  1.23  1.23;
%               2.45  0.5   0.5;
%               2.01  0.7   0.7;
%               2.25  0.6   0.6;
%               1.9   0.8   0.8];
% 
% 
% %
% % tot_ves = max(writeconn(:));
% % tot_term = length(terminal);
% 
% % Write all this information to file.
% dlmwrite('connectivity.txt',writeconn,'\t');
% dlmwrite('terminal_vessels.txt',terminal,'\t');
% dlmwrite('dimensions.txt',dimensions,'\t');

%%
writeconn  = dlmread('data_fluids/Sample0/connectivity.txt');
terminal   = dlmread('data_fluids/Sample0/terminal_vessels.txt');
dimensions = dlmread(append('data_fluids/Sample', num2str(sample_number), '/dimensions.txt'));
 
tot_ves = max(writeconn(:))+1;
tot_term = length(terminal);

Qdat    = dlmread('Qdat_SimVes8192.dat');
scale   = 1;
Qdat    = Qdat.*scale;
dlmwrite('Qin.dat',Qdat,'delimiter','\t')

%% Parameter Values
% MJC
f1   = 2.5e+6;%5e+6;
f2   = -15; %10
f3   = 6.4e+4;%8e+4;%8e+4;

fs1 = 10*f1;
fs2 = f2;
fs3 = 8e5;
% fs1  = 1e5;%f1;%3*f1;%2.5e+7;%;%5e+6;%f1;
% fs2  = -14;%-15;%-20;%f2;
% fs3  = 5e4;%f3;%3*f3;%8e6;%1e+6;%f3;

% Structured tree values (used to determine radius of each vessel, 
% r = r_root * alpha^m * beta^n)
% alpha, beta < 1, alpha + beta > 1
alpha = 0.88;%0.84;
beta  = 0.697;%0.7;  

% Length to radius ratio
lrr = 15;

% Minimum radius (cm)
r_min = 0.008;%0.001;


% Number points evaluated along each large vessel, 
% smallest length > 1/num_pts 
num_pts = 6;%24;

% File ID (useful for running the C++/Fortran code in parallel)
ID = 1;

tot_sten = 0;
tot_screen = 0;
%% Put the parameters into a vector and run the model
par_nom = [f1, f2, f3, fs1, fs2, fs3,...
           alpha, beta, lrr, r_min,...
           tot_ves, tot_term, tot_sten, tot_screen, num_pts, ID, sample_number];
       
param_str = mat2str(par_nom);

% Run the model
% NOTE: Windows users need 'sor06.exe', Mac/Linux users need ./sor06
tic
out = unix(sprintf('./sor06 %s -w',param_str(2:end-1)));
toc
if out == 1
    disp 'there is a model output'
else
    disp 'there is no model output'
end

%% Load all the model results
% NOTE: Results are stored such that the first 1:N entries are the PROXIMAL
% large vessel predictions, N+1:2N are the
% MIDPOINT predictions in the same vessels, and 2N+1:3N are the DISTAL
% predictions.

name = sprintf('output_%d.2d',ID);

data = load(name);

%% Extract data of interest

% t - time 
% x - length (cm)
% p - pressure
% q - flow
% a - area
% c - wave speed
[~,x,p,q,a,c] = gnuplot(data); % extract data

nv = tot_ves; % total number of vessels

ntp = size(p,1); % number of time points in the time series

pressure_in = p(:,1:nv); % pressure at the inlet of vessels
pressure_mid = p(:,nv+1:2*nv); % pressure at the midpoint of vessels
pressure_out = p(:,2*nv+1:3*nv); % pressure at the outlet of vessels

flow_in = q(:,1:nv); % flow at the inlet of vessels
flow_mid = q(:,nv+1:2*nv); % flow at the midpoint of vessels
flow_out = q(:,2*nv+1:3*nv); % flow at the outlet of vessels

area_in = a(:,1:nv); % area at the inlet of vessels
area_mid = a(:,nv+1:2*nv); % area at the midpoint of vessels
area_out = a(:,2*nv+1:3*nv); % area at the outlet of vessels

% Time vector
T = 0.85; % Period
t = linspace(0,T,ntp);

% Velocity - u = q/a
u = q./a;

u_in = flow_in./area_in;
u_mid = flow_mid./area_mid;
u_out = flow_out./area_out;

omega = 2*pi/T; % rad
mu    = 0.032;  % g/cm/s - viscosity 
rho   = 1.055;  % g/ml - density 
nu    = mu/rho; % cm^2/s - kinematic viscosity
delta = sqrt(nu/omega); % cm

shear_stress_in = mu*u_in./delta;
shear_stress_mid = mu*u_mid./delta;
shear_stress_out = mu*u_in./delta;
%% Do some plotting

% % Plot inlet pressure and flow
% figure(1);clf(1);
% plot(t,pressure_mid(:, 1:3),'LineWidth',3);
% ylabel('Pressure (mmHg)');
% xlabel('Time (s)');
% grid on; set(gca,'FontSize',20);
% 
% % Plot inlet pressure and flow
% figure(2);clf(2)
% plot(t,flow_mid(:, 1:3),'LineWidth',3);
% ylabel('Flow (mL/s)')
% xlabel('Time (s)');
% grid on; set(gca,'FontSize',20);
% 
% %%
% figure;
% plot(t,shear_stress_mid,'LineWidth',3)
% ylabel('Shear Stress (g/cm \cdot s^2)')
% xlabel('Time (s)');
% grid on; set(gca,'fontsize',20)

%% Save data

save(append('flow_in_output/flow_in', num2str(sample_number), '.mat'), 'flow_in')
save(append('pressure_in_output/pressure_in', num2str(sample_number), '.mat'), 'pressure_in')
save('flow_in_output/time_step.mat', "t")
save('pressure_in_output/time_step.mat', "t")

save(append('flow_mid_output/flow_mid', num2str(sample_number), '.mat'), 'flow_mid')
save(append('pressure_mid_output/pressure_mid', num2str(sample_number), '.mat'), 'pressure_mid')
save('flow_mid_output/time_step.mat', "t")
save('pressure_mid_output/time_step.mat', "t")

save(append('flow_out_output/flow_out', num2str(sample_number), '.mat'), 'flow_out')
save(append('pressure_out_output/pressure_out', num2str(sample_number), '.mat'), 'pressure_out')
save('flow_out_output/time_step.mat', "t")
save('pressure_out_output/time_step.mat', "t")

end