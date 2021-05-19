function WIA(t,p,q,a,r,stiff,plot_ID)
% Define constants
rho = 1.055;
n = length(p);
p = p.*1333.22;
u = q./a;
dp = diff(p);
du = diff(u);
dt = t(2)-t(1);
tnew = t(1:end-1);
% NOTE: W/m^2 is 0.001 g/s^3
Wconv = 0.001;

% Fudge factor: C was not nondimensionalized correctly
c = wave_speed(a,r,stiff);

rho_cW = rho.*c(1:end-1);

dp_fwd = 0.5.*(dp+rho_cW.*du);
dp_bwd = 0.5.*(dp-rho_cW.*du);

P_fwd = fwd_bwd(dp_fwd,p);
P_bwd = fwd_bwd(dp_bwd,p);

du_fwd = 0.5.*(du+(dp./rho_cW));
du_bwd = 0.5.*(du-(dp./rho_cW));

U_fwd = fwd_bwd(du_fwd,u);
U_bwd = fwd_bwd(du_bwd,u);

WI_fwd = Wconv.*(dp_fwd.*du_fwd)./dt^2;
WI_bwd = Wconv.*(dp_bwd.*du_bwd)./dt^2;

% figure(plot_ID);clf;
% subplot(2,2,1); hold on;
% plot(tnew,dp_fwd./1333.22,'r');
% plot(tnew,dp_bwd./1333.22,'b');
% ylim([-0.05 0.1]);
% set(gca,'FontSize',24);
% xlim([0 tnew(end)]);
% 
% subplot(2,2,2); hold on;
% plot(tnew,P_fwd./1333.22,':r');
% plot(tnew,P_bwd./1333.22,':b');
% plot(t,p./1333.22,'k');
% set(gca,'FontSize',24);
% xlim([0 tnew(end)]);
% 
% subplot(2,2,3); hold on;
% plot(tnew,du_fwd,'r');
% plot(tnew,du_bwd,'b');
% ylim([-0.4 0.6]);
% set(gca,'FontSize',24);
% xlim([0 tnew(end)]);
% 
% subplot(2,2,4); hold on;
% plot(tnew,U_fwd,':r');
% plot(tnew,U_bwd,':b');
% plot(t,u,'k');
% set(gca,'FontSize',24);
% xlim([0 tnew(end)]);

% figure; plot(dp_fwd,'r'); hold on; plot(dp_bwd,'b');


% cw = tnew(1:id);
% dw = tnew((id+1):n-1);
% figure(plot_ID+1);hold on;
% 
% h = area(cw',WI_fwd(1:id)); h.FaceColor = 'r';
% h = area(dw',WI_fwd(id+1:end)); h.FaceColor = 'c';
% 
% h = area(cw',WI_bwd(1:id)); h.FaceColor = 'b';
% h = area(dw',WI_bwd(id+1:end)); h.FaceColor = 'm';
% ylim([-3.5e4 8.2e4]);
% xlim([0 tnew(end)]);
% set(gca,'FontSize',24);

% disp([sum(P_bwd(1:id))./sum(P_fwd(1:id)) max(P_bwd)./max(P_fwd)]);

%% Try plotting differently
ids = 1:length(dp_fwd);
bool_fwd = dp_fwd>0;
bool_bwd = dp_bwd>0;
fwd_comp = ids(bool_fwd);
fwd_exp  = ids(~bool_fwd);
bwd_comp = ids(bool_bwd);
bwd_exp  = ids(~bool_bwd);
%%
figure(plot_ID+1);hold on;
h = area(tnew(fwd_comp),WI_fwd(fwd_comp)); 
h.FaceColor = [0.1 0.1 0.1]; %h.LineWidth = 2;
h = area(tnew(fwd_exp),WI_fwd(fwd_exp)); 
h.FaceColor = [0.4 0.4 0.4]; %h.LineWidth = 2;
h = area(tnew(bwd_comp),WI_bwd(bwd_comp)); 
h.FaceColor = [0.6 0.6 0.6]; h.LineWidth = 3; h.LineStyle = ':';
h = area(tnew(bwd_exp),WI_bwd(bwd_exp)); 
h.FaceColor = [0.9 0.9 0.9]; h.LineWidth = 3; h.LineStyle = ':';
plot(tnew,0.*tnew,'k','LineWidth',3);

% figure(plot_ID+1);hold on;
% plot(fwd_comp,WI_fwd(fwd_comp),'or','LineWidth',2,'MarkerSize',10)
% plot(fwd_exp,WI_fwd(fwd_exp),'oc','LineWidth',2,'MarkerSize',10)
% plot(bwd_comp,WI_bwd(bwd_comp),'ob','LineWidth',2,'MarkerSize',10)
% plot(bwd_exp,WI_bwd(bwd_exp),'om','LineWidth',2,'MarkerSize',10)
ylim([-3.5e4 8.2e4]);
xlim([0 tnew(end)]);
set(gca,'FontSize',24);

disp([sum(P_bwd(fwd_comp))./sum(P_fwd(fwd_comp)) max(P_bwd)./max(P_fwd)]);
end

function c = wave_speed(a,r,fr)
a0 = pi.*r.^2;
c =  sqrt(0.5.*fr.*sqrt(a/a0));
% c =  0.5*fr(r)*sqrt(a/a0); % use if fr is a function
end

function out = fwd_bwd(dX,X)
n = length(dX);
out = zeros(n,1);
out(1) = dX(1);
for i=2:n
    out(i) = out(i-1)+dX(i);
end
out = out + X(1);
end

