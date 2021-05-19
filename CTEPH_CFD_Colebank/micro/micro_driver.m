% Run the one-sided structured tree model and predict pressure and flow in
% the small arteries
%
% MJC 7/10/2020
clear; clc; %close all;
%
% load ../scenario_b_FIXED.mat;
% load ../scenario_c_V4.mat;
% load ../scenario_d_V4.mat;
load ../scenario_e_FIXED.mat
tot_ves = size(p,2)/3;
% If dealing with stenosis scenario, add this line

%% Want to define an area function
stiff = @(r0,f1,f2,f3) f1.*exp(f2.*r0)+f3;
A = @(p,r0,f1,f2,f3) (pi.*r0.^2).*(1+(3./4./stiff(r0,f1,f2,f3)).*p).^2;
%%
p_h = cell(length(terminal),2);
q_h = cell(length(terminal),2);
a_h = cell(length(terminal),2);
D_h = cell(length(terminal),2);

counter = 1;

%% For terminals that are blocked, we need to increase z0
% Only apply for scenarios c,d
load('../term_conn_sten.mat','block','screen');
path = get_vessel_path(conn,block,screen);


% Decide which vessels to use
block_ids   = [];
unblock_ids = 1:length(terminal);
for path_id = 1:size(path,1)
    for i=1:length(path{path_id})
        if any(path{path_id}(i)==terminal)
            block_ids(end+1) = find(path{path_id}(i)==terminal);
        end
    end
end
for id = terminal
% id = terminal(end-1);
%%
p0 = p(:,id+tot_ves*2);
q0 = q(:,id+tot_ves*2);
tmpts = length(p0);
r_root = dim_mat(id,3);
Period = 0.85;

f1   = pars(1);%5e+6;
f2   = pars(2); %10
f3   = pars(3);
fs1  = pars(4);%5e+7;%f1;
fs2  = pars(5);%-20;%f2;
fs3  = pars(6);%5e+6*1;%f3;
alpha = pars(7); %Alpha
beta = pars(8); %Beta
lrr  =  pars(9);
r_min   = pars(10);
if any(id==terminal(block_ids))
    trm_rst = 1e8; %% CHANGE FOR SCENARIO e
else
    trm_rst = pars(11);
end
pars_ST = [Period trm_rst fs1 fs2 fs3 alpha beta lrr r_min r_root]';
dlmwrite('parameters.dat',pars_ST,'delimiter','\t','precision',16);
dlmwrite('pterm.dat',p0,'precision',16);

run_tree

%%
% Determine the number of plots generated
alpha_r = r_root.*alpha.^(0:100);
beta_r = r_root.*beta.^(0:100);
n = find(alpha_r<r_min,1)-1;% alpha generations
m = find(beta_r<r_min,1)-1; % beta generations

Palpha = zeros(tmpts,n);
Qalpha = zeros(tmpts,n);
Time = zeros(tmpts,n);
Pbeta = zeros(tmpts,m);
Qbeta = zeros(tmpts,m);


%%
ra_d = [];
rb_d = [];

for i = 1:n
    s = num2str(i-1);
    s = strcat('p',s,'_alpha.2d');
    data = load(s);
    Palpha(:,i) = data(:,1);
    Qalpha(:,i) = data(:,2);%*10/T; % Dimensionalize correctly!
%     Time_a(:,i) = t;
    ra_d = [ra_d r_root*(alpha^(i-1))];
end

for i = 1:m
    s = num2str(i-1);
    s = strcat('p',s,'_beta.2d');
    data = load(s);
    Pbeta(:,i) = data(:,1);
    Qbeta(:,i) = data(:,2);  % Dimensionalize correctly!
%     Time_b(:,i) = t;
    rb_d = [rb_d r_root*(beta^(i-1))];
end


% Predict area from pressure

pa_ND = Palpha.*1333.22;
pb_ND = Pbeta.*1333.22;

A_alpha = A(pa_ND,ra_d,f1,f2,f3);
A_beta = A(pb_ND,rb_d,f1,f2,f3);
D_alpha = (A_alpha-(pi.*ra_d.^2))./(pi.*ra_d.^2);
D_beta = (A_beta-(pi.*rb_d.^2))./(pi.*rb_d.^2);

% figure(10); hold on; plot(linspace(0,1,size(Palpha,2)),mean(Palpha)); 
% figure(20); hold on;  plot(linspace(0,1,size(Pbeta,2)),mean(Pbeta)); 
% figure(30); hold on; plot(linspace(0,1,size(Qalpha,2)),mean(Qalpha)); 
% figure(40); hold on;  plot(linspace(0,1,size(Qbeta,2)),mean(Qbeta)); 
% 
% figure(id); 
% subplot(2,2,1); plot(Palpha);hold on; plot(p0,'k')
% subplot(2,2,2); plot(Pbeta);hold on; plot(p0,'k')
% subplot(2,2,3); plot(Qalpha);hold on; plot(q0,'k')
% subplot(2,2,4); plot(Qbeta);hold on; plot(q0,'k')

p_h{counter,1} = Palpha;
p_h{counter,2} = Pbeta;
q_h{counter,1} = Qalpha;
q_h{counter,2} = Qbeta;
a_h{counter,1} = A_alpha;
a_h{counter,2} = A_beta;
D_h{counter,1} = D_alpha;
D_h{counter,2} = D_beta;
counter = counter+1;

mu = 0.032;
T = 0.85;
rho = 1.055;
spars = [alpha,beta,T,mu,rho,fs1,fs2,fs3,r_root,r_min];
% out = get_shear_micro(Palpha,Qalpha,Pbeta,Qbeta,spars);
end
save('se_ST_trial','p_h','q_h','a_h','D_h');