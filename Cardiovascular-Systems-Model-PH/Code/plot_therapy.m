%% Call CV model and plot optimial prediction
clear; clc;
load('PostAll.mat','Xopt');
%%
pat = 2;
load('Patient_2/P2Shift_02_P2_ALL.mat');
data = data{1};
pars_OG = history{1}.parsOPT;
history = DriverBasic_postOPT(data,pars_OG,5,pat);
[p,V,q] = CVmodel(pars_OG,data);

%%
pars = pars_OG;
pars(2) = log(exp(pars(2)).*0.75);
pars(11) = log(exp(pars(11)).*1.25);
history = DriverBasic_postOPT(data,pars,1,pat);
[p1,V1,q1] = CVmodel(pars,data);

%%
pars = pars_OG;
pars(2) = log(exp(pars(2)).*0.5);
pars(11) = log(exp(pars(11)).*1.5);
history = DriverBasic_postOPT(data,pars,3,pat);
[p2,V2,q2] = CVmodel(pars,data);

%%
pars = pars_OG;
pars(2) = log(0.0646);
pars(11) = log(4.1815);
history = DriverBasic_postOPT(data,pars,2,pat);
[p3,V3,q3] = CVmodel(pars,data);

%% Now make a plot
figure(1); clf;
subplot(1,3,1); hold on;
plot(p(:,5),'r','LineWidth',3)
plot(p1(:,5),'c','LineWidth',3)
plot(p2(:,5),'m','LineWidth',3)
plot(p3(:,5),'b','LineWidth',3)
set(gca,'FontSize',28);
axis tight; grid on; xticklabels({});
ylabel('Pressure (mmHg)');
title('RA');

subplot(1,3,2); hold on;
plot(p(:,6),'r','LineWidth',3)
plot(p1(:,6),'c','LineWidth',3)
plot(p2(:,6),'m','LineWidth',3)
plot(p3(:,6),'b','LineWidth',3)
set(gca,'FontSize',28);
axis tight; grid on; xticklabels({});
title('RV');

subplot(1,3,3); hold on;
plot(p(:,7),'r','LineWidth',3)
plot(p1(:,7),'c','LineWidth',3)
plot(p2(:,7),'m','LineWidth',3)
plot(p3(:,7),'b','LineWidth',3)
set(gca,'FontSize',28); xticklabels({});
axis tight; grid on;
title('PA');

%% Look at 9 panel plot

figure(100); %clf;
subplot(3,3,1); hold on;
plot(p(:,1),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('pla (mmHg)');

subplot(3,3,2);hold on;
plot(p(:,2),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('plv (mmHg)');

subplot(3,3,3);hold on;
plot(p(:,3),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('psa (mmHg)');

subplot(3,3,4);hold on;
plot(p(:,4),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('psv (mmHg)');

subplot(3,3,5); hold on;
plot(p(:,5),'r','LineWidth',3.2);
 set(gca,'FontSize',24);
 ylabel('pra (mmHg)');


subplot(3,3,6); hold on;
plot(p(:,6),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('prv (mmHg)');



subplot(3,3,7);hold on;
plot(p(:,7),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('ppa (mmHg)');
xlabel('time (s)');
    
tdc = linspace(0,1,100);
subplot(3,3,8);hold on;
plot(tdc,mean(ppv)*ones(size(tdc)),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('ppv (mmHg)');
xlabel('time (s)')
ylim([mean(ppv)*.9 , mean(ppv)*1.1])

CO = trapz(q(:,3));
subplot(3,3,9);hold on;
plot(tdc,CO*ones(size(tdc)),'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('CO (mL/sec)');
xlabel('time (s)')


