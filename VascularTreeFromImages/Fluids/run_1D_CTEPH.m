%=============================================================*
%                                                             *
% run_1D_ST.m                                                 *
% Version: 1.1 (created on 2 Apr. 2020)                       *
% AUTHORS: M.J. Colebank, M.U. Qureshi,  M.S. Olufsen         *
%          Department of Mathematics                          *
%          North Carolina State University, Raleigh, USA      *
% DATE UPDATED: Apr. 2, 2020.                                 *
%                                                             *
% DESCRIPTION: This script creates the inteface with the C++  *
% code for the 1D fluid dynamics model. The script creates    * 
% and runs the executable by passing selected model           *
% parameters and plots the hemodynamic waveforms              *
%                                                             *
%=============================================================*



clear; close all; clc;
! make
%%
f = @(q,r) q(1).*exp(q(2).*r) + q(3);
f1   = 2.5e+6;%5e+6;
f2   = -15; %10
f3   = 1e+5;%8e+4;%8e+4;
fs1  = f1;
fs2  = f2;
fs3  = f3;
Z0  = 0;
num_pts = 12;
ID = 1;
%% Scale stiffness

alpha = 0.88; %Alpha
beta = 0.66; %Beta
lrr  =  50; %20
r_min   = 0.001; % Diameter of 15 microns at terminal radius

%% Add a stenosis 
% Right now, ring like lesions are added at junction points in the parent
% branch; ask MJC for different code if you want lesions in the middle of a
% branch.

% Ring lesions
sten_val = 0.50;% This controls what percent occlusion you'd like
max_LL = 0.25;  % This controls what portion of the vessel (%) the lesion occupies
% Web lesions (don't recommend using unless interested in CTEPH)
permeability = 0; % Degree of web like lesion severity (0.01 --> 1e-4, the smaller the more severe)
screen_L     = 0.0; % This controls what portion of the vessel (%) the lesion occupies

% This dictates WHERE the lesion is in the network. Make these two empty
% arrays if no lesions are necessary
block       = [2];
screen      = [];

%% Upscale flow
scale = 1.1;
qdat = load('Qdat_SimVes8192.dat');
qdat = qdat.*scale;
dlmwrite('Qin.dat',qdat);

% Converging full network
writeconn  	 = [0 1 2 0
                1 3 4 0
                2 5 6 0];
% %             
terminal = [3 4 5 6];


dimensions = [3.58  1.27  1.27;
              6.24  1.19  1.19;
              5.58  1.23  1.23;
              2.45  0.5   0.5;
              2.01  0.7   0.7;
              2.25  0.6   0.6;
              1.9   0.8   0.8];
%% Stenosis part
sten_factor = ones(1,length(block)).*sten_val;
sten_length = dimensions(block,1).*max_LL.*sten_val;
screen_factor = ones(1,length(screen)).*permeability;
screen_length = dimensions(screen,1).*screen_L;

sten_L = 1./num_pts;
if ~isempty(block) || ~isempty(screen)
    if ~isempty(screen)
        [conn, terminal, dimensions,block,screen] = make_stenosis_screen(writeconn,terminal,dimensions,block,screen,screen_L,num_pts);
        writeconn = conn;
        ids = writeconn>0;
        writeconn = writeconn-ids;
    else
        writeconn = make_stenosis_2(writeconn,block-1);
    end
end
tot_ves = max(writeconn(:))+1;
tot_term = length(terminal);
sten_write_tofile = [block'-1 sten_factor' sten_length];
screen_write_tofile = [screen'-1 screen_factor' screen_length];
%% Write everything to file
% load scenario_e_V4.mat
dlmwrite('connectivity.txt',writeconn,'\t');
dlmwrite('terminal_vessels.txt',terminal,'\t');
dlmwrite('Dimensions.txt',round(dimensions,2),'\t')
dlmwrite('Sten.txt',sten_write_tofile,'\t');
dlmwrite('Screen.txt',screen_write_tofile,'\t');
%% Nominal
tot_sten   = length(block);
tot_screen = length(screen);
%% Put the parameters into a vector and run the model
par_nom = [f1, f2, f3, fs1, fs2, fs3,...
           alpha, beta, lrr, r_min,...
           tot_ves, tot_term, tot_sten, tot_screen, num_pts, ID];
       
param_str = mat2str(par_nom);

% Run the model
% NOTE: Windows users need 'sor06.exe', Mac/Linux users need ./sor06
tic
out = unix(sprintf('./sor06 %s',param_str(2:end-1)));
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

% p - pressure
% q - flow
[~,~,p,q,~,~] = gnuplot(data); % extract data

nv = tot_ves; % total no of vessels

ntp = size(p,1); % no of time points in the flow or pressure time series

pressure_all = (p(:,nv+1:2*nv)); % middle prediction all vessels
flow_all = (q(:,nv+1:2*nv)); % middle prediction all vessels



%% Do some plotting
% Time vector
t = linspace(0,0.85,ntp);

% Plot midpoint in the vessels
figure(1);clf(1);
plot(t,p,'LineWidth',3);
ylabel('Pressure (mmHg)');
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);

figure(2);clf(2)
plot(t,q,'LineWidth',3);
ylabel('Flow (mL/s)')
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);