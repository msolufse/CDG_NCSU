% Code for manuscript "Parameter inference in a computational model of 
% hemodynamics in pulmonary hypertension" by Colebank, Colunga, et al.
% This driver file will call each patient's optimized parameters and plot
% the model predictions compared to the static and dynamics data.

%
clear; clc; close all;
%%
clear; clc;
PATIENT   = 1;  %% Change for each patient (1-5), or use 6 for normotensive
TREATMENT = 7; %% Change for each treatment (1-10)
[p,V,q]    = Driver_postopt(PATIENT,2);
[pT,VT,qT] = Driver_treatment(PATIENT,TREATMENT); % Returns pressure, volume, flow

%% Now compare pressure predictions
figure(100); clf;

subplot(2,2,1); hold on;
plot(p(:,2),'r','LineWidth',3);
plot(pT(:,2),'--b','LineWidth',3);
set(gca,'FontSize',20);
ylabel('p_{lv} (mmHg)');
hold off;

subplot(2,2,2); hold on;
plot(p(:,3),'r','LineWidth',3);
plot(pT(:,3),'--b','LineWidth',3);
set(gca,'FontSize',20);
ylabel('p_{sa} (mmHg)');
hold off;


subplot(2,2,3); hold on;
plot(p(:,6),'r','LineWidth',3);
plot(pT(:,6),'--b','LineWidth',3);
set(gca,'FontSize',20);
ylabel('p_{rv} (mmHg)');
hold off;

subplot(2,2,4); hold on;
plot(p(:,7),'r','LineWidth',3);
plot(pT(:,7),'--b','LineWidth',3);
set(gca,'FontSize',20);
ylabel('p_{pa} (mmHg)');
hold off;