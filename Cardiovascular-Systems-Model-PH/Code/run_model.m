% Code for manuscript "Parameter inference in a computational model of 
% hemodynamics in pulmonary hypertension" by Colebank, Colunga, et al.
% This driver file will call each patient's optimized parameters and plot
% the model predictions compared to the static and dynamics data.

%

%% Call CV model and plot optimial prediction
clear; clc;
PATIENT = 1;  %% Change for each patient (1-5), or use 6 for normotensive
RESIDUAL = 2; %% Change for each residual (1-2)
[p,V,q] = Driver_postopt(PATIENT,RESIDUAL); % Returns pressure, volume, flow