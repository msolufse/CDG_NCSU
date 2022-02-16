close all;
clc;
clear variables;
%% Decide which vessel to plot
load control_LargeVesselResults.mat

which_mouse = 2;
load(strcat('C',num2str(which_mouse)));
connectivity = conn{which_mouse};
g = geometry{which_mouse};
p = pressures{which_mouse};
q = flows{which_mouse};
parent = connectivity(:,1)+1;
d1 = connectivity(:,2)+1;
d2 = connectivity(:,3)+1;
term = [];

figure(99); clf; hold on;
% First plot the whole network
arcs = arcs_p1C2;
for i=1:size(arcs,2)
    plot3(-1.*arcs{1,i}(2:end,1),1.*arcs{1,i}(2:end,2),-1.*arcs{1,i}(2:end,3),'o','Color',[0.9 0.9 0.9]);
end
% Now plot the principal pathway
for i=1:length(parent)
    plot3(-g{parent(i),1}(:,1),g{parent(i),1}(:,2),-g{parent(i),1}(:,3),'ok');
    if ~any(d1(i)==parent)
        if length(g{d1(i),1})>1
            term(end+1) = d1(i);
            plot3(-g{d1(i),1}(:,1),g{d1(i),1}(:,2),-g{d1(i),1}(:,3),'or');
        end
    end
    if ~any(d2(i)==parent)
        if length(g{d2(i),1})>1
            term(end+1) = d2(i);
            plot3(-g{d2(i),1}(:,1),g{d2(i),1}(:,2),-g{d2(i),1}(:,3),'or');
        end
    end
end
%%
num_ves = length(g)-1;
r_mat = [];
for i=term
    r_mat(end+1) = g{i,2};
end
id = find(r_mat < 0.3);
% Find pathway TO this terminal vessel
% for pick_id = 30:length(id)
% pick_id = 31;
min_diff = 0;
for pick_id=31%1:length(id)
    plot3(-g{term(id(pick_id)),1}(:,1),g{term(id(pick_id)),1}(:,2),-g{term(id(pick_id)),1}(:,3),'oc');
    path = [term(id(pick_id))];
    new_path = path;
    r_principal = [];
    while ~isempty(new_path)
        id_path_d1 = find(d1==path(end));
        id_path_d2 = find(d2==path(end));
        new_path = [id_path_d1 id_path_d2];
        if ~isempty(new_path)
            path(end+1) = new_path;
            r_principal(end+1) = g{new_path,2};
        end
    end
    r_principal = fliplr(r_principal);
    if min(r_principal)>=0.028
        disp([pick_id length(r_principal)])
    end
figure(8); clf; hold on;
        plot(r_principal,'--o');
        disp([pick_id length(r_principal) r_mat(id(pick_id))])
        
        figure(10); clf; hold on;
        plot(mean(p(:,path)),'-o');
end
% pick_id = pick_id_keep;
p_principal = p(:,num_ves+path);
q_principal = q(:,num_ves+path);

% Added by MJC, 3/11/2020
% Need to find the radius of the principal pathway for plotting purposes
r_principal = [];
for i=path
    r_principal(end+1) = g{i,2};
end
%
r_root = r_mat(id(pick_id))
p_in = p(:,num_ves + term(id(pick_id)));
q_in = q(:,num_ves + term(id(pick_id)));
figure(222); plot(p_in);
% end
t = linspace(0,0.11,1024);
dlmwrite('p_mouse.dat',[t;p_in']','delimiter','\t');
%%
run_tree
%% Parameters below should match parameters.dat
rho = 1.055;
% Radius dependent length
lrr    = @(r) (13.39.*exp(-0.007708.*r.*10.^4)).*r;
% lrr = 15.75;

T = 0.11;
t = 0:T/1023:T;
mu = 0.049;

aa  = 0.8837;
bb = 0.6666;
r_min = 0.005;

% Determine the number of plots generated
alpha_r = r_root.*aa.^(1:100);
beta_r = r_root.*bb.^(1:100);
n = find(alpha_r<r_min,1);% alpha generations
m = find(beta_r<r_min,1); % beta generations

Palpha = zeros(1024,n);
Qalpha = zeros(1024,n);
Time = zeros(1024,n);
Pbeta = zeros(1024,m);
Qbeta = zeros(1024,m);


%%
ra_d = [];
rb_d = [];

for i = 1:n
    s = num2str(i-1);
    s = strcat('p',s,'_alpha.2d');
    data = load(s);
    Palpha(:,i) = data(:,1);
    Qalpha(:,i) = data(:,2)*10/T; % Dimensionalize correctly!
    Time_a(:,i) = t;
    ra_d = [ra_d r_root*(aa^(i-1))];
end

for i = 1:m
    s = num2str(i-1);
    s = strcat('p',s,'_beta.2d');
    data = load(s);
    Pbeta(:,i) = data(:,1);
    Qbeta(:,i) = data(:,2)*10/T;  % Dimensionalize correctly!
    Time_b(:,i) = t;
    rb_d = [rb_d r_root*(bb^(i-1))];
end


La = lrr(ra_d);
Lb = lrr(rb_d); % vessel's length

Lac = cumsum(La);
Lbc = cumsum(Lb); % cumulative pathe length


%% Start plots
%% Pressure and Flow
color_space_alpha = [1:-1./(n+2):0; (0:1./(n+2):1).*0; (0:1./(n+2):1).*0]';
color_space_beta  = [(0:1./(m+2):1).*0; (0:1./(m+2):1).*0; 1:-1./(m+2):0]';
color_space_prin  = [0:1./(length(path)+5):1; 0:1./(length(path)+5):1; 0:1./(length(path)+5):1]';
% Plot Principal pathway in grey shades
figure (1); clf; hold on;
for i=1:length(path)
   plot(t,p_principal(:,i),'LineWidth',3,'Color',color_space_prin(i,:))
end
% Plot Alpha branch
for i=1:n
plot(t,Palpha(:,i),'LineWidth',3,'Color',color_space_alpha(i,:))
end

for i=1:m
plot(t,Pbeta(:,i),'LineWidth',3,'Color',color_space_beta(i,:))
end
set(gca,'fontsize',30,'xlim',[0 T]);
ylabel 'P (mmHg)'
xlabel 't (s)'
ylim([0 30])
print -depsc Pall.eps

% Just alpha and Beta branches
figure (100); clf; hold on;
% Plot Alpha branch
for i=1:n
plot(t,Palpha(:,i),'LineWidth',3,'Color',color_space_alpha(i,:))
end
set(gca,'fontsize',30,'xlim',[0 T]);%, 'ylim', [-5 30]);
ylabel 'P (mmHg)'
xlabel 't (s)'
ylim([0 8])
% title '\alpha branch'
print -depsc Palpha_Z0.eps




figure (200);  clf; hold on;
for i=1:m
plot(t,Pbeta(:,i),'LineWidth',3,'Color',color_space_beta(i,:))
end
set(gca,'fontsize',30,'xlim',[0 T]);%, 'ylim', [-5 30]);
ylabel 'P (mmHg)'
xlabel 't (s)'
ylim([0 8])
% title '\beta branch'
print -depsc Pbeta_Z0.eps
%%
%% Pressure and Flow
% Plot Principal pathway in grey shades
figure (10); clf; hold on;
for i=1:length(path)
   plot(t,q_principal(:,i),'LineWidth',3,'Color',color_space_prin(i,:))
end
% Plot Alpha branch
for i=1:n
plot(t,Qalpha(:,i),'LineWidth',3,'Color',color_space_alpha(i,:))
end

for i=1:m
plot(t,Qbeta(:,i),'LineWidth',3,'Color',color_space_beta(i,:))
end
set(gca,'fontsize',30,'xlim',[0 T]);
ylabel 'Q (ml/s)'
xlabel 't (s)'
ylim([0 0.45])
print -depsc Qall.eps




figure (11); clf;hold on;
for i=1:n
plot(t,Qalpha(:,i),'LineWidth',3,'Color',color_space_alpha(i,:))
end
set(gca,'fontsize',30,'xlim',[0 T]);%, 'ylim', [-0.2 1.8]);
ylabel 'Q (ml/s)'
xlabel 't (s)'
ylim([0 7e-3])
% title '\alpha branch'
print -depsc Qalpha_Z0.eps

figure (12); clf; hold on;
for i=1:m
plot(t,Qbeta(:,i),'LineWidth',3,'Color',color_space_beta(i,:))
end
set(gca,'fontsize',30,'xlim',[0 T]);%, 'ylim', [-0.2 1.8]);
ylabel 'Q (ml/s)'
xlabel 't (s)'
ylim([0 7e-3])
% title '\beta branch'
print -depsc Qbeta_Z0.eps


%% Plot steady components against radius
% MJC 3/11/2020, add mean components from principal pathway
figure (4); clf;
% semilogy(r_principal,mean(p_principal),'--o','LineWidth',2.5,'Color',[0.7 0.7 0.7]);
hold on;
semilogy(ra_d,mean(Palpha),'-or','linewidth',2.5);
semilogy(rb_d,mean(Pbeta),'--ob','linewidth',2.5);
%plot(linspace(r_principal(1),r_principal(end),length(r_principal)),mean(p_principal),'--o','LineWidth',2.5,'Color',[0.7 0.7 0.7]);
set(gca,'Fontsize',30,'Xdir','Reverse');
% legend('\alpha branch', '\beta branch','Location','SouthWest')
xlabel 'radius (cm)'
ylabel 'Mean P (mmHg)'
ylim([2.5 4.5])
grid on;
print -depsc Pmean_Z0.eps

figure (5); clf;
% semilogy([r_principal(1) r_principal],[mean(Qalpha(:,1)) mean(q_principal)],'--o','LineWidth',2.5,'Color',[0.7 0.7 0.7]);
hold on;
semilogy(ra_d,mean(Qalpha),'-or','linewidth',2.5); 
semilogy(rb_d,mean(Qbeta),'--ob','linewidth',2.5);
set(gca,'Fontsize',30,'Xdir','Reverse');
% legend('\alpha branch', '\beta branch','Location','SouthWest')
xlabel 'radius (cm)'
ylabel 'Mean Q (ml/s)'
ylim([0 4.5e-3])
grid on;    
print -depsc Qmean_Z0.eps

% figure (6);
% plot(ra_d,mean(Tau_alpha),'-or','linewidth',3); hold on;
% plot(rb_d,mean(Tau_beta),'--ob','linewidth',3);
% set(gca,'Fontsize',20,'Xdir','Reverse');
% legend('\alpha branch', '\beta branch','Location','NorthWest')
% xlabel 'radius (cm)'
% ylabel '\tau_{avg}'
% print -depsc Taumean_Z0.eps

%% Despite huge difference in the number of generations max and min radii of both branches is same
% so ploting against the radius vector without markers may give a false impression about the behavior
% To better see the behavior, here we plot them against cumulative path
% length of the branch (beta is much shorter than the alpha branch)
figure (40);
plot(Lac,mean(Palpha),'-or','linewidth',3); hold on;
plot(Lbc,mean(Pbeta),'--ob','linewidth',3);
set(gca,'Fontsize',20);
legend('\alpha branch', '\beta branch','Location','best')
xlabel 'Path length (cm)'
ylabel 'Pmean (mmHg)'
print -depsc Pmean_pathlength_Z0.eps

figure (50);
plot(Lac,mean(Qalpha),'-or','linewidth',3); hold on;
plot(Lbc,mean(Qbeta),'--ob','linewidth',3);
set(gca,'Fontsize',20);
legend('\alpha branch', '\beta branch','Location','best')
xlabel 'Path length (cm)'
ylabel 'Qmean (ml)'
print -depsc Qmean_pathlength_Z0.eps

% figure (60);
% plot(Lac,mean(Tau_alpha),'-or','linewidth',3); hold on;
% plot(Lbc,mean(Tau_beta),'--ob','linewidth',3);
% set(gca,'Fontsize',20);
% legend('\alpha branch', '\beta branch','Location','best')
% xlabel 'Path length (cm)'
% ylabel '\tau_{avg}'
% print -depsc Taumean_pathlength_Z0.eps

%% Calculate the reynolds number
q_principal = fliplr(q_principal);
r_principal = fliplr(r_principal);
f1   = 1e+4;%9e+3;
f2   = -50;
f3   = 3.6e+4;%3.9e+4;
fs1  = 5e+4;
fs2  = -60;
fs3  = 2.0e+4;
fstiffLA = @(r) (4./3).*f1.*exp(f2.*r) + f3;
fstiffST = @(r) (4./3).*fs1.*exp(fs2.*r) + fs3;
afuncLA  = @(p,r0,A0) A0./(1-p*1333.22./fstiffLA(r0)).^2; 
afuncST  = @(p,r0,A0) A0./(1-p*1333.22./fstiffST(r0)).^2; 
diameter_a = [r_principal ra_d].*2;%.*0.1;
dmicro_a = diameter_a.*10^4;

diameter_b = [r_principal rb_d].*2;%.*0.1;
dmicro_b = diameter_b.*10^4;

% Alpha side
eta_D = 3.2+6.*exp(-0.085.*dmicro_a)-2.44.*exp(-0.06.*(dmicro_a.^0.645));
C_exp = (0.8+exp(-0.075.*dmicro_a)).*(1.0./(1.0+10.^(-11.0) .* dmicro_a.^12.0) - 1.0)...
        + (1.0./(1.0+10.^(-11.0) .* dmicro_a.^12.0));
Hct = 0.81;
Hct_term = ((1.0 - Hct).^C_exp - 1.0)/((1.0-0.45).^C_exp - 1.0);
mu_Da  = (1.0 + (eta_D - 1.0).*Hct_term.*...
       (dmicro_a./(dmicro_a - 1.1)).^2.0).*(dmicro_a./(dmicro_a - 1.1)).^2.0;
mu1a = mu.*mu_Da./3.2;

% Beta side
eta_D = 3.2+6.*exp(-0.085.*dmicro_b)-2.44.*exp(-0.06.*(dmicro_b.^0.645));
C_exp = (0.8+exp(-0.075.*dmicro_b)).*(1.0./(1.0+10.^(-11.0) .* dmicro_b.^12.0) - 1.0)...
        + (1.0./(1.0+10.^(-11.0) .* dmicro_b.^12.0));
Hct = 0.81;
Hct_term = ((1.0 - Hct).^C_exp - 1.0)/((1.0-0.45).^C_exp - 1.0);
mu_Db  = (1.0 + (eta_D - 1.0).*Hct_term.*...
       (dmicro_b./(dmicro_b - 1.1)).^2.0).*(dmicro_b./(dmicro_b - 1.1)).^2.0;
mu1b = mu.*mu_Db./3.2;

A_alpha = afuncST(Palpha,ra_d,pi.*ra_d.^2);
A_beta = afuncST(Pbeta,rb_d,pi.*rb_d.^2);

A_princ = afuncLA(p_principal,r_principal,pi.*r_principal.^2);
A_pa = [A_princ A_alpha];
A_pb = [A_princ A_beta];


u_alpha = Qalpha./A_alpha;
u_beta  = Qbeta./A_beta;
u_princ = q_principal./A_princ;
u_pa  = [u_princ u_alpha];
u_pb  = [u_princ u_beta];

Ra = rho.*u_pa.* diameter_a./ mu1a; 
Ra(:,1:length(r_principal)) = rho.*u_princ.* diameter_a(1:length(r_principal))./ mu;

Rb = rho.*u_pb.* diameter_b./ mu1b; 
Rb(:,1:length(r_principal)) = rho.*u_princ.* diameter_b(1:length(r_principal))./ mu;
