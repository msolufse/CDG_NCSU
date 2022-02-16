%% Run the ST model with parameters
% MJC 6/1/2019'
clear; %close all; clc;
! make clean -f Makefile
! make -f Makefile
! chmod +x sor06

%% Define parameters
f = @(q,r) q(1).*exp(q(2).*r) + q(3);
num_generations = 13;
% PARAMETERS
f1   = 1e+4;
f2   = -50;
f3   = 3.6e+4;

fs1  = 5e+4;
fs2  = -60;
fs3  = 2.0e+4;

rm   = 0.005;


pressures = {};
flows     = {};
geometry  = {};
ST        = {};
%% Set the flow profile
PQ_ID = 6;
qdat = load(strcat('qC',num2str(PQ_ID),'_8192.dat'));
pdat = load(strcat('pC',num2str(PQ_ID),'_1024.dat'));
qdat = qdat;
dlmwrite('QIN.dat',qdat);
BIFSCALE = [12   10 7.4;...
    13.8 13 8.8; ...
    11.4 10.6 9];
%% Nominal
for i=3%1:3
    filename = strcat('PP_C',num2str(i),'_p1.mat');
    load(filename);
    [geo3d,tot_ves,tot_term] = create_data_fluids(num_generations, ves_conn, term, dim_mat,net3D,BIFSCALE(i,:));%,r_std);
    alpha = 0.8837;
    beta  = 0.6666;
    lrr1  = 13.39;
    lrr2  = -0.007708;
    pars = [f1 f2 f3 fs1 fs2 fs3 alpha beta lrr1 lrr2 rm tot_ves tot_term];
    num_pars = length(pars);
    
    %% Plot the exponential stiffness
    D = load('Dimensions.txt');
    figure; hold on;
    rspace = linspace(0,max(D(:,2)),1000);
    plot(rspace,f([f1 f2 f3],rspace),'r');
    plot(D(:,2),f([f1 f2 f3],D(:,2)),'k*','LineWidth',2);
    %% Call the model
    sc = []; %Only if you want to do optimization and use scaling
    [p,q,a,c] = STmodel_V2(pars,sc);
    %% Plot the results
    t = linspace(0,0.11,1024);
    if ~isempty(p)
        for j=1:3
            figure(10 + j); hold on;
            plot(t,p(:,j),'LineWidth',3);
            plot(t,p(:,tot_ves+j),'LineWidth',3);
            figure(1000+j); hold on;
            plot(t,q(:,j),'LineWidth',3);
            plot(t,q(:,tot_ves+j),'LineWidth',3);
        end
    end
    pressures{end+1} = p;
    flows{end+1}     = q;
    geometry{end+1}  = geo3d;
    ST{end+1}        = plot_ST_resistance(tot_term);
    pressure_drop = max(p(:,1))-max(p(:,tot_ves+1))
end

%% Load data
meanQ = [];
PQ_ID = 6
qdat = load(strcat('qC',num2str(PQ_ID),'_8192.dat'));
pdat = load(strcat('pC',num2str(PQ_ID),'_1024.dat'));
figure(11); hold on;
plot(t,pdat,'k','LineWidth',3);

tlong = linspace(0,0.11,8193);
figure(1001); hold on;
plot(tlong,qdat,'k','LineWidth',3);


