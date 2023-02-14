% Code for manuscript "Parameter inference in a computational model of 
% hemodynamics in pulmonary hypertension" by Colunga, Colebank, REU, and Olufsen
% This driver file will call each patient's optimized parameters and plot
% the model predictions compared to the static and dynamics data.

% There are TWO adjustable variables for this driver script
% PATIENT - Change for each PH patient (1-9) or for normotensive (10)
% RESIDUAL - Change for either (1) static residual or (2) static + dynamic

%% Call CV model and plot optimial prediction
clear; clc;
PATIENT = 6;  %% Change for each patient (1-5), or use 6 for normotensive
RESIDUAL = 2; %% Change for each residual (1-2)
Call_Model(PATIENT,RESIDUAL); % Returns pressure, volume, flow