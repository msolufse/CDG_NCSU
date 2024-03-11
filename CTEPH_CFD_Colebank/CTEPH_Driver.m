%=============================================================*
%                                                             *
% CTEPH_Driver.m                                              *
% Version: 2.0 (created on 21 July. 2021)                     *
% AUTHORS: M.J. Colebank, M.U. Qureshi,  M.S. Olufsen         *
%          Department of Mathematics                          *
%          North Carolina State University, Raleigh, USA      *
% DATE UPDATED: March 5, 2024.                                *
%                                                             *
% DESCRIPTION: This script creates the inteface with the C++  *
% code for the 1D fluid dynamics model. The script creates    *
% and runs the executable by passing selected model           *
% parameters and plots the hemodynamic waveforms. The         *
% driver reproduces the conditions and results from the       *
% manuscript "A multiscale model of vascular function in      *
% chronic thromboembolic pulmonary hypertension" by           *
% Colebank et al., 2021                                       *
%=============================================================*
%%
clear; close all; clc;
! make clean -f Makefile
! make -f Makefile
! chmod +x sor06

%% NOTE: 
% Zc flag can be set to 0 if you do not want to recalculate the structured
% tree parameters, since this takes time. Only set Zc_flag=0 if you are
% keeping the structured tree parameters FIXED.
Zc_flag = 1;
if Zc_flag==1 % Otherwise, erase all the STnetwork files and recalculate the impedances
    ! make veryclean
    ! make
end
%%
% Define parameters and scenarios
% Scenario index 1-5 corresponds to normotension, a-d in th manuscript
scenario = 5;

if scenario == 1 || scenario == 2 % Control OR Lesions only (norm and case a)
    f1   = 2.5e+6;
    f2   = -15;
    f3   = 6.4e+4;
    fs1  = 2.5e7;
    fs2  = -15;
    fs3  = 8e5;
    Z0  = 1e5;
elseif scenario == 3 || scenario == 4 % Hypertensive with no PVR OR with PVR increase (case b and c)
    f1   = 2.5e+6;
    f2   = -20;
    f3   = 8e+4;
    fs1  = 5e7;
    fs2  = -20;
    fs3  = 1.6e6;
    if scenario == 3
        Z0  = 1e6; % CTEPH - no PVR increase - case b
    else
        Z0  = 1e7;% CTEPH - PVR inrease - case c
    end
elseif scenario==5 % Severe CTEPH with large vessel dilation (case d)
    f1   = 1e+7;
    f2   = -10;
    f3   = 2e5;
    fs1  = 5e7;
    fs2  = -20;
    fs3  = 1.6e6;
    Z0  = 1e8;
else
    error('Scenario not found. Exiting')

end

% Set up the locations for lesions
if scenario==1
    block      = [];
    screen     = [];
    sten_val   = 0.0;
    max_LL     = 0.0;
    screen_val = 0.0;
    screen_L   = 0.0;  
else
    block       = [14 20 33 35 51 65 93 95 102];
    screen      = [17 23 29 36 50 56 63 75 89 100 107];
    sten_val   = 0.90;
    max_LL     = 0.25;
    screen_val = 0.009;
    screen_L   = 0.50;  
end
%% Other parameters in the system (does not change)
alpha = 0.88; 
beta  = 0.68; 
lrr  =  0; % Note, the LRR is fixed in the root_imp.f90 file
rm   = 0.00075; % Diameter of 15 microns at terminal radius


%% Generate the network file
generations = 14;
num_pts  = 8; % Spatial resolution in PDEs is 1/num_pts
[~,~,conn,dim_mat,terminal,geo3D] = create_data_fluids2(generations,0,0,num_pts);
write_conn = max(conn-1,0);

%% Generate the stenosis files
sten_factor = ones(1,length(block)).*sten_val;
sten_length = dim_mat(block,1).*max_LL.*sten_val;
screen_factor = ones(1,length(screen)).*screen_val;
screen_length = dim_mat(screen,1).*screen_L;

if ~isempty(block) || ~isempty(screen)
    if ~isempty(screen)
        [conn, terminal, dim_mat,block,screen] = make_stenosis_screen(conn,terminal,dim_mat,block,screen,screen_L,num_pts);
        write_conn = conn;
        ids = write_conn>0;
        write_conn = write_conn-ids;
    else
        write_conn = make_stenosis_2(write_conn,block-1);
    end
end
tot_ves = max(conn(:));     tot_term = length(terminal);
sten_write_tofile = [block'-1 sten_factor' sten_length];
screen_write_tofile = [screen'-1 screen_factor' screen_length];
%% Only for cases b-d
% Note: this is nested conditionals, since case d has all the features of
% cases a, b, and c, e.g.
vessel_stiff_vec = [];
vessel_stiff = [];
if scenario > 2
    % If we want to increase resistance in small vessels, we need to know all
    % the vessels downstream
    path = get_vessel_path(conn,block,screen);
    % Consolodate all these into a vector of vessels to constrict
    lesions = [block screen];
    for i=1:length(lesions)
        n_temp = length(path{i});
        temp = path{i};
        vessel_stiff(end+1:end+n_temp) = temp;
        id = find(vessel_stiff == lesions(i));
        vessel_stiff(id)=[];
        % DO NOT NARROW VESSELS IN A WEB
        if any(lesions(i)==screen)
            id = find(vessel_stiff == lesions(i)+1);
            vessel_stiff(id)=[];
        end
    end
    vessel_stiff = unique(vessel_stiff);

    % Narrow branches distal to lesion by 40%
    narrow = 0.6;
    r_old = dim_mat(:,2);
    dim_mat(vessel_stiff,2:3) = dim_mat(vessel_stiff,2:3).*narrow;

    if scenario>3 % Narrow all branches less than 5mm2 area
        R_mm = dim_mat(:,2)*10; % Convert to mm
        CSA  = pi.*R_mm.^2;     % Calculate area
        not_sten = 1:tot_ves; not_sten(vessel_stiff) = []; % Only narrow branches not affected by case b above
        ids_narrow = not_sten(CSA(not_sten)<5);
        dim_mat(ids_narrow,2:3) = dim_mat(ids_narrow,2:3).*0.6;

        if scenario==5 % Dilate proximal arteries as well as the above
            % Widening
            ids_widen = not_sten(CSA(not_sten)>10);
            dim_mat(ids_widen,2:3) = dim_mat(ids_widen,2:3).*1.1;
        end
    end
    n_stiff = length(vessel_stiff);
    vessel_stiff_vec = [vessel_stiff'-1 rm.*ones(n_stiff,1) alpha*ones(n_stiff,1) beta*ones(n_stiff,1)];
end

%% Write everything to file
dlmwrite('connectivity.txt',write_conn,'\t');
dlmwrite('terminal_vessels.txt',terminal-1,'\t');
dlmwrite('Dimensions.txt',round(dim_mat,3),'\t')
dlmwrite('Sten.txt',round(sten_write_tofile,3),'\t');
dlmwrite('Screen.txt',screen_write_tofile,'\t');
dlmwrite('stiff.txt',vessel_stiff_vec,'\t');
%% Define parameter vector
tot_stiff  = length(vessel_stiff);
tot_sten   = length(block);
tot_screen = length(screen);
pars = [f1 f2 f3 fs1 fs2 fs3 alpha beta lrr rm Z0 tot_stiff tot_ves ...
    tot_term tot_sten tot_screen num_pts Zc_flag];
num_pars = length(pars);

%% Call the model (change to './sor06' if using Mac or Linux)
clc;
str_pars = mat2str(pars);
call = unix(sprintf('sor06.exe %s',str_pars(2:end-1)));
%% Plot the results
names = {};
load('pu_ALL.2d'); % Stores the solutions at the proximal, midpoint, and distal location
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
    end
end

