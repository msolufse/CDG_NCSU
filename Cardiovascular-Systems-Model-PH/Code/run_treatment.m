% Code for manuscript "Parameter inference in a computational model of 
% hemodynamics in pulmonary hypertension" by Colebank, Colunga, et al.
% This driver file will call each patient's optimized parameters and plot
% the model predictions compared to the static and dynamics data.

%
clear; clc; close all;
%%
clear; clc;
PATIENT   = 1;  %% Change for each patient (1-5), or use 6 for normotensive
TREATMENT = 1; %% Change for each treatment (1-10)
Driver_treatment(PATIENT,TREATMENT);