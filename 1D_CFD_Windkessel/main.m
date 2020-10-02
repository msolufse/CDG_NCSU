% Driver file for a single vessel model of pulse wave propagation.
% Model usines a simple linear elastic wall model and a three element Windkessl outflow
% boundary condition. The system is driven by imposing an inflow profile.
% For this example parameters were taken from a recent study by Flores et al.
% See Table 2 in Annals of Biomedical Engineering, Vol. 44, No. 10. pp.
% 3047?3068 (2016) DOI: 10.1007/s10439-016-1625-3

clc;
close all;
clear all;

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create an Inflow waveform for pulse wave simlation in a single vessel

Qmin = 0;    % minumum flow rate during diastole (ml/s)
Qmax = 500;  % maximum flow rate during systole (ml/s)


t0 = 0.0; % initial time (s)
T = 0.85; % Length of the cardiac cycle (s)

Tm  = 0.15; % Peak systole time (s)
Td  = 0.4;  % End systole time (s)

tt = Td-Tm; % time difference between peak and end of systole

[Qin] = InFlow(Qmin,Qmax,t0,Tm,tt,T); % Generates a desired inflow profile and saves it to a file Qin_*.dat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify Geometric and heomodynamic parameter (Human Thoracic Aorta) 

Length = 24.14; % (cm) Length of Human Thoracic Aorta
Rin = 1.2;   % (cm) Inlet radius
Rout = 1.2;  % (cm) Outlet radius. (Reduce rout marginally to introduce tapering)
r0 = mean([Rin Rout]); % (cm) A reference radius to claculate stiffness.
                       % Set to either rin or rout for straight vessels

CO = mean(Qin); % Cardiac output (ml/s)
Pmean = 97.0;   % Desired mean pressure (mmHg)

Inductance = 0.0; % Inductance of the 4-element Windkessel Model. 
                  % Set zero for three element Windkessel Model (mmHg.s^2/ml)                
tau = 1.313;      % Time constant of exponentional diastolic decay
alpha = 0.14;     % Ratio of R1 to RT in the Windkessel Model.
                  % Set alpha = 1 to reduce to a 2-element Windkessel Model

h = 0.12; %  Wal thickness (cm)
E = 3800; %  Constant Young's Modulus (mmHg). Increase to make the vessel stiffer.

HB = 12;     % Number of cardiac cycles ran to acheive steady sate
cycles = 1;  % Integer number of cycles you want plot in graphs
id = 1;      % Change it to save output files corresponding to different parameter combinition
nbrves = 1;  % Number of vessels. !! DONOT CHANGE. THE CODE IS HARD CODED FOR SINGLE VESSEL ONLY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Following function computes non-dimensiol arterial stiffenss and Windkessel parameters

[Ehr0, R1, R2, CT, L] = ND_Par(E, h, r0, CO, Pmean, alpha, tau, Inductance);

% Pass parameters to C++ excuitable 'sor06'

unix(sprintf('./sor06 %f %f %f %f %f %f %f %f %d %d %d %d',...
                Length, Rin, Rout, Ehr0, R1, R2,...
                CT, L, nbrves, HB, cycles, id));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting Simulated data

data = load ('pu1_1.2d');

[t,x,p,q,A,C] = gnuplot(data);

P_p = p(:,1);  Q_p = q(:,1);   A_p = A(:,1); C_p = C(:,1);
P_m = p(:,floor((end)/2));  Q_m = q(:,floor((end)/2));   A_m = A(:,floor((end)/2)); C_m = C(:,floor((end)/2));
P_d = p(:,end); Q_d = q(:,end); A_d = A(:,end); C_d = C(:,end);
t = t(:,1)-(HB-cycles)*0.85;

Qmean_sim = [mean(Q_p) mean(Q_d)]
Pmean_sim = [mean(P_p) mean(P_d)]

figure(2);
plot(t,P_p,t,P_m,t,P_d,'linewidth',2)
set(gca, 'fontsize',20);grid on;
xlabel 'Time (s)';axis tight;
ylabel 'Pressure (mmHg)';
legend ('Proximal', 'mid', 'distal','Location','Best')

figure(3);
plot(t,Q_p,t,Q_m,t,Q_d,'linewidth',2)
set(gca, 'fontsize',20);grid on;
xlabel 'Time (s)';axis tight;
ylabel 'Flow (ml/s)';
legend ('Proximal', 'mid', 'distal','Location','Best')

figure(4);
plot(Q_p,P_p,Q_m,P_m,Q_d,P_d,'linewidth',2)
set(gca, 'fontsize',20);grid on;
xlabel 'Flow (ml/s)';
ylabel 'Pressure (mmHg)';
legend ('Proximal', 'mid', 'distal','Location','Best')

figure(5);
plot(P_p,A_p,P_m,A_m,P_d,A_d,'linewidth',2)
set(gca, 'fontsize',20);grid on;
xlabel 'Pressure (mmHg)';
ylabel 'Area (cm^2)';
legend ('Proximal', 'mid', 'distal','Location','Best')
