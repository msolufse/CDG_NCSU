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



clear; %close all; clc;
! make clean -f Makefile
! make -f Makefile
! chmod +x sor06
% clc;
% Define parameters
f = @(q,r) q(1).*exp(q(2).*r) + q(3);

% % Control
f1   = 2.5e+6;%5e+6;
f2   = -15; %10
f3   = 8e+4;%8e+4;%8e+4;
fs1  = f1*10;%5e+6;%f1;
fs2  = f2;%-20;%f2;
fs3  = f3*10;%1e+6;%f3;
Z0  = 1e5;
f3 = f3.*0.8;

% Hypertensive
% f1   = 2.5e+6;%5e+6;
% f2   = -14; %10
% f3   = 8e+4;%8e+4;
% fs1  = f1*20;%5e+6;%f1;
% fs2  = f2;%-20;%f2;
% fs3  = f3*20;%1e+6;%f3;


% Z0  = 1e5;% CTEPH - no PVR increase - scenario c
% Z0  = 1e7;% CTEPH - PVR inrease - scenario d
% Z0  = 1e8; %CTEPH -PVR increase - scenario e

% MJC new addition for scenario e
% f1   = 1e+7;%2.5e+6;%5e+6;
% f2   = -10;%-14; %10
% f3   = f3*2.5;%1.8;%8e+4;

%% Scale stiffness


par1 = 0.88; %Alpha
par2 = 0.68; %Beta
lrr  =  20; %20
rm   = 0.00075; % Diameter of 15 microns at terminal radius

sten_val = 0.50;%0.80;
max_LL = 0.25;
permeability = 0.009;%5e-1;
screen_L     = 0.50;
% %% Upscale flow
scale = 1.6;
qdat = load('Qdat_SimVes8192.dat');
qdat = qdat.*scale;
qdat = filter_signal(qdat,10,8192*2+1);
dlmwrite('QIN.dat',qdat');

%% Use a bigger network, MJC 3/25/2020
generations = 14;
Zc_flag = 1;%1;
if Zc_flag==1
    ! make veryclean
    ! make
end
num_pts  = 8;%8;%10;
[tot_ves,tot_term,conn,dim_mat,terminal,geo3D] = create_data_fluids2(generations,0,0,num_pts);


ftemp1 = f([f1 f2 f3],dim_mat(:,2));
rspace = linspace(rm,min(dim_mat(:,2)),100);
ftemp2 = f([fs1 fs2 fs3],rspace);
figure(1);clf; semilogy(dim_mat(:,2),ftemp1,'*')
hold on; semilogy(rspace,ftemp2,'r*');


figure(2); clf; hold on;
for i=1:length(geo3D)
   plot3(geo3D{i}(:,1),geo3D{i}(:,2),geo3D{i}(:,3),'LineWidth',3) 
end
tot_ves = max(conn(:));
tot_term = length(terminal);
write_conn = max(conn-1,0);


% Calculate the terminal impedance
% Pterm = 5;
% r_term = dim_mat(terminal,2);
% Z0  =calculate_Z0(Pterm,qdat,tot_term,r_term,rm,par1,par2);
% Z0   = 8e+10;


%% Add a stenosis to one of the terminal vessels
block       = [14 20 33 35 51 65 93 95 102];
screen      = [17 23 29 36 50 56 63 75 89 100 107];


sten_factor = ones(1,length(block)).*sten_val;
sten_length = dim_mat(block,1).*max_LL.*sten_val;
screen_factor = ones(1,length(screen)).*permeability;
screen_length = dim_mat(screen,1).*screen_L;

sten_L = 1./num_pts;
if ~isempty(block) || ~isempty(screen)
% [conn, terminal, dim_mat] = make_stenosis(conn,terminal,dim_mat,sten_loc,sten_L,num_pts);
% write_conn = make_stenosis_2(write_conn,block-1); % Subtract one to match connectivity file
    if ~isempty(screen)
        [conn, terminal, dim_mat,block,screen] = make_stenosis_screen(conn,terminal,dim_mat,block,screen,screen_L,num_pts);
        write_conn = conn;
        ids = write_conn>0;
        write_conn = write_conn-ids;
    else
        write_conn = make_stenosis_2(write_conn,block-1);
    end
end
tot_ves = max(conn(:));
tot_term = length(terminal);
sten_write_tofile = [block'-1 sten_factor' sten_length];
screen_write_tofile = [screen'-1 screen_factor' screen_length];
%%
% If we want to increase resistance in small vessels, we need to know all
% the vessels downstream
path = get_vessel_path(conn,block,screen);
% Consolodate all these into a vector of vessels to constrict
vessel_stiff = [];
lesions = [block screen];
for i=1:length(lesions)
    n_temp = length(path{i});
    temp = path{i};
    vessel_stiff(end+1:end+n_temp) = temp;
    id = find(vessel_stiff == lesions(i));
    vessel_stiff(id)=[];
    % DO NOT NARROW VESSELS IN A WEB
%     if any(lesions(i)==screen)
%         id = find(vessel_stiff == lesions(i)+1);
%         vessel_stiff(id)=[];
%     end
end
vessel_stiff = unique(vessel_stiff);
%% Scenario C:
% Add this line to introduce vessel narrowing to vessels distal a blockage
% narrow = 0.6;
% r_old = dim_mat(:,2);
% dim_mat(vessel_stiff,2:3) = dim_mat(vessel_stiff,2:3).*narrow;
% 
% % % Add this as well to change rmin to match
% rmin_new = zeros(length(vessel_stiff),1);
% alpha_new = zeros(length(vessel_stiff),1);
% beta_new = zeros(length(vessel_stiff),1);
% for i=1:length(vessel_stiff)
%     if any(vessel_stiff(i)==terminal)
%         rmin_new(i) = rm*1;%*narrow;
%         alpha_new(i) = par1*1;%narrow;
%         beta_new(i) = par2*1;%narrow;
%     end
% end
% vessel_stiff = [vessel_stiff'-1 rmin_new alpha_new beta_new];

%% Add this line to introduce vessel widening for A>10mm
% % and narrowing in other vessels
%% Use this to increase vessel narrowing for cross-sec area < 5mm
% R_mm = dim_mat(:,2)*10;
% CSA  = pi.*R_mm.^2;
% not_sten = 1:tot_ves;
% not_sten(vessel_stiff(:,1)+1) = [];


% Narrowing
% ids_narrow = not_sten(CSA(not_sten)<5);
% dim_mat(ids_narrow,2:3) = dim_mat(ids_narrow,2:3).*0.6;
% 
% % Widening
% ids_widen = not_sten(CSA(not_sten)>10);
% dim_mat(ids_widen,2:3) = dim_mat(ids_widen,2:3).*1.1;


%%% DOES NOT WORK
% % vessel_stiff(end+1:end+length(ids_widen),:) = [ids_widen'-1 zeros(length(ids_widen),1) zeros(length(ids_widen),1)  zeros(length(ids_widen),1) ];
% % [who,where] = sort(vessel_stiff(:,1),'ascend');
% % vessel_stiff = vessel_stiff(where,:);

%% Add this line to introduce vessel widening (where there is no stenosis)
% not_sten = 1:tot_ves;
% not_sten(vessel_stiff(:,1)+1) = [];
% dim_mat(not_sten,2:3) = dim_mat(not_sten,2:3).*1.2;
% f3 = 1.5*f3;
%% Use this to increase vessel narrowing EVERYWHERE
% not_sten = 1:tot_ves;
% not_sten(vessel_stiff(:,1)+1) = [];
% dim_mat(not_sten,2:3) = dim_mat(not_sten,2:3).*0.80;
%% Use this to increase vessel narrowing for diameter<500microns
% R_cuttoff = 0.05/2;
% not_sten = 1:tot_ves;
% not_sten(vessel_stiff(:,1)+1) = [];
% ids = not_sten(dim_mat(not_sten,2)<=R_cuttoff);
% dim_mat(ids,2:3) = dim_mat(ids,2:3).*0.8;

% %% Use this to increase vessel narrowing for cross-sec area < 5mm
% R_mm = dim_mat(:,2)*10;
% CSA  = pi.*R_mm.^2;
% not_sten = 1:tot_ves;
% not_sten(vessel_stiff(:,1)+1) = [];
% ids = not_sten(CSA(not_sten)<5);
% dim_mat(ids,2:3) = dim_mat(ids,2:3).*0.8;
%% Use this to induce vascular recruitment
% par1 = par1.*0.9; par2 = par2.*1.3;
% par1 = 0.9; par2 = 0.82;
% rm = rm*0.5;
%% Use this to increase vessel narrowing in terminal vessels only
% not_sten = 1:tot_ves;
% not_sten(vessel_stiff(:,1)+1) = [];
% not_sten = not_sten(any(not_sten==terminal'));
% dim_mat(not_sten,2:3) = dim_mat(not_sten,2:3).*0.80;
% % Or do nothing
vessel_stiff = [];
%% Write everything to file
dlmwrite('connectivity.txt',write_conn,'\t');
dlmwrite('terminal_vessels.txt',terminal-1,'\t');
dlmwrite('Dimensions.txt',round(dim_mat,2),'\t')
dlmwrite('Sten.txt',sten_write_tofile,'\t');
dlmwrite('Screen.txt',screen_write_tofile,'\t');
dlmwrite('stiff.txt',vessel_stiff,'\t');
%% Nominal
tot_stiff  = size(vessel_stiff,1);
tot_sten   = length(block);
tot_screen = length(screen);
pars = [f1 f2 f3 fs1 fs2 fs3 par1 par2 lrr rm Z0 tot_stiff tot_ves ...
        tot_term tot_sten tot_screen num_pts Zc_flag];
num_pars = length(pars);

%% Call the model
clc;
str_pars = mat2str(pars);
call = unix(sprintf('./sor06 %s',str_pars(2:end-1)));
%% Plot the results
names = {};
load('pu_ALL.2d');
[~,~,p,q,a,c] = gnuplot(pu_ALL);
t = linspace(0,0.85,size(p,1));
if call==0
    for j=1:3%tot_ves
        figure(10 + j); hold on;
        plot(t,p(:,j),'LineWidth',3);
        plot(t,p(:,tot_ves+j),'--','LineWidth',3);
        legend('Prox 1','Prox 2','Dist 1','Dist 2')
        figure(1000+j); hold on;
        plot(t,q(:,j),'LineWidth',3);
        plot(t,q(:,tot_ves+j),'LineWidth',3);
        names{end+1} = strcat('V',num2str(j));
        
%         figure(2000+j); hold on;
%         plot(t,a(:,j),'LineWidth',3);
%         plot(t,a(:,tot_ves+j),'LineWidth',3);
%         names{end+1} = strcat('V',num2str(j));
    end
end


% ST = plot_ST_resistance(tot_term);
pressure_drop = max(p(:,1))-max(p(:,tot_ves+1))
