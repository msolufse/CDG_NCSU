%% Perform wave intensity analysis
clear; clc; close all;

sa = load('scenario_a_V2.mat');
sb = load('scenario_b_FIXED.mat');
sc = load('scenario_c_FIXED.mat');
sd = load('scenario_d_FIXED.mat');
sf = load('scenario_e_FIXED.mat');
% sf = load('scenario_f_5mm_narrow_10mm_dilate_Zc_ALL.mat');
% sf = load('scenario_f2_Zc.mat');
t = linspace(0,0.85,1024);
%%
tot_ves = size(sa.p,2)/3;
f1   = sa.pars(1);%5e+6;
f2   = sa.pars(2); %10
f3   = sa.pars(3);%8e+4;
r = sa.dim_mat(1:3,2);
fa = f1.*exp(f2.*r) + f3;
% WIA(t,sa.p(:,1),sa.q(:,1),sa.a(:,1),r,fa,100)
WIA(t,sa.p(:,1+tot_ves),sa.q(:,1+tot_ves),sa.a(:,1+tot_ves),r(1),fa(1),110)
WIA(t,sa.p(:,2+tot_ves),sa.q(:,2+tot_ves),sa.a(:,2+tot_ves),r(2),fa(2),120)
WIA(t,sa.p(:,3+tot_ves),sa.q(:,3+tot_ves),sa.a(:,3+tot_ves),r(3),fa(3),130)

%%
tot_ves = size(sb.p,2)/3;
f1   = sb.pars(1);%5e+6;
f2   = sb.pars(2); %10
f3   = sb.pars(3);%8e+4;
r = sb.dim_mat(1:3,2);
fb = f1.*exp(f2.*r) + f3;

% WIA(t,sb.p(:,1),sb.q(:,1),sb.a(:,1),r,fa,200)
WIA(t,sb.p(:,1+tot_ves),sb.q(:,1+tot_ves),sb.a(:,1+tot_ves),r(1),fb(1),210)
WIA(t,sb.p(:,2+tot_ves),sb.q(:,2+tot_ves),sb.a(:,2+tot_ves),r(2),fb(2),220)
WIA(t,sb.p(:,3+tot_ves),sb.q(:,3+tot_ves),sb.a(:,3+tot_ves),r(3),fb(3),230)


%%
tot_ves = size(sc.p,2)/3;
f1   = sc.pars(1);%5e+6;
f2   = sc.pars(2); %10
f3   = sc.pars(3);%8e+4;
r    = sc.dim_mat(1:3,2);
fc   = f1.*exp(f2.*r) + f3;

WIA(t,sc.p(:,1+tot_ves),sc.q(:,1+tot_ves),sc.a(:,1+tot_ves),r(1),fc(1),310)
WIA(t,sc.p(:,2+tot_ves),sc.q(:,2+tot_ves),sc.a(:,2+tot_ves),r(2),fc(2),320)
WIA(t,sc.p(:,3+tot_ves),sc.q(:,3+tot_ves),sc.a(:,3+tot_ves),r(3),fc(3),330)

%%
%%
tot_ves = size(sd.p,2)/3;
f1   = sd.pars(1);%5e+6;
f2   = sd.pars(2); %10
f3   = sd.pars(3);%8e+4;
r    = sd.dim_mat(1:3,2);
fd   = f1.*exp(f2.*r) + f3;

WIA(t,sd.p(:,1+tot_ves),sd.q(:,1+tot_ves),sd.a(:,1+tot_ves),r(1),fd(1),410)
WIA(t,sd.p(:,2+tot_ves),sd.q(:,2+tot_ves),sd.a(:,2+tot_ves),r(2),fd(2),420)
WIA(t,sd.p(:,3+tot_ves),sd.q(:,3+tot_ves),sd.a(:,3+tot_ves),r(3),fd(3),430)



%%
% tot_ves = size(se.p,2)/2;
% r = se.dim_mat(1,2);
% fe = f1.*exp(f2.*r) + f3;
% 
% % WIA(t,se.p(:,1),se.q(:,1),se.a(:,1),r,fe,500)
% WIA(t,se.p(:,1+tot_ves),se.q(:,1+tot_ves),se.a(:,1+tot_ves),r,fe,510)



%%
tot_ves = size(sf.p,2)/3;
f1   = sf.pars(1);%5e+6;
f2   = sf.pars(2); %10
f3   = sf.pars(3);%8e+4;
r    = sf.dim_mat(1:3,2);
ff   = f1.*exp(f2.*r) + f3;

WIA(t,sf.p(:,1+tot_ves),sf.q(:,1+tot_ves),sf.a(:,1+tot_ves),r(1),ff(1),510)
WIA(t,sf.p(:,2+tot_ves),sf.q(:,2+tot_ves),sf.a(:,2+tot_ves),r(2),ff(2),520)
WIA(t,sf.p(:,3+tot_ves),sf.q(:,3+tot_ves),sf.a(:,3+tot_ves),r(3),ff(3),530)


%% Network level analysis

% f1   = 2.5e+6;%5e+6;
% f2   = -15; %10
% f3   = 6.4e+4;%8e+4;
% stiff_f = @(rin) f1.*exp(f2.*rin) + f3;
% network_WIA(t,sa.p(:,1:floor(end/2)),sa.q(:,1:floor(end/2)),sa.a(:,1:floor(end/2)),sa.dim_mat(:,2),stiff_f,1000)
% 
% 
% f1 = sf.pars(1);
% f2 = sf.pars(2);
% f3 = sf.pars(3);
% stiff_f = @(rin) f1.*exp(f2.*rin) + f3;
% 
% network_WIA(t,sf.p(:,1:floor(end/2)),sf.q(:,1:floor(end/2)),sf.a(:,1:floor(end/2)),sf.dim_mat(:,2),stiff_f,2000)

%%
figure(111); print('WIA/sa_1','-dpdf');
figure(121); print('WIA/sa_2','-dpdf');
figure(131); print('WIA/sa_3','-dpdf');
figure(211); print('WIA/sb_1','-dpdf');
figure(221); print('WIA/sb_2','-dpdf');
figure(231); print('WIA/sb_3','-dpdf');
figure(311); print('WIA/sc_1','-dpdf');
figure(321); print('WIA/sc_2','-dpdf');
figure(331); print('WIA/sc_3','-dpdf');
figure(411); print('WIA/sd_1','-dpdf');
figure(421); print('WIA/sd_2','-dpdf');
figure(431); print('WIA/sd_3','-dpdf');
figure(511); print('WIA/se_1','-dpdf');
figure(521); print('WIA/se_2','-dpdf');
figure(531); print('WIA/se_3','-dpdf');