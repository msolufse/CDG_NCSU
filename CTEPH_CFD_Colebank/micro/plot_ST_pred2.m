%% Load in structured tree predictions and plot
clear; clc; close all;
load control_ST.mat
load sb_ST.mat
p_hb = p_h;
q_hb = q_h;
a_hb = a_h;

% load sc_ST.mat
% p_hc = p_h;
% q_hc = q_h;
% a_hc = a_h;
load se_ST.mat
load term_IDs.mat
load ../term_conn_sten.mat

%%
t = linspace(0,0.85,1024);
% Loop through blocked vessels
for i=[10,15, 25]%length(block_ids)
figure(i+100); clf;

subplot(2,3,1); 
plot(t,p_c{block_ids(i),1}(:,1:2:end),'-r','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,2);
plot(t,p_hb{block_ids(i),1}(:,1:2:end),'-','Color',[0.82 0.52 0.63],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,3);
plot(t,p_h{block_ids(i),1}(:,1:2:end),'-','Color',[1.0 0.5 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,4); hold on;
plot(t,p_c{block_ids(i),2},'-b','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,5);
plot(t,p_hb{block_ids(i),2},'-','Color',[0.1 1.0 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,6);
plot(t,p_h{block_ids(i),2},'-','Color',[0.0 0.9 0.9],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);
print(strcat('../../Results/ST/p_block_',num2str(i),'.pdf'),'-dpdf')

%% Flow

figure(i+300); clf;
subplot(2,3,1); 
plot(t,q_c{block_ids(i),1}(:,1:2:end),'-r','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,2);
plot(t,q_hb{block_ids(i),1}(:,1:2:end),'-','Color',[0.82 0.52 0.63],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,3);
plot(t,q_h{block_ids(i),1}(:,1:2:end),'-','Color',[1.0 0.5 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,4); hold on;
plot(t,q_c{block_ids(i),2},'-b','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,5);
plot(t,q_hb{block_ids(i),2},'-','Color',[0.1 1.0 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,6);
plot(t,q_h{block_ids(i),2},'-','Color',[0.0 0.9 0.9],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);
print(strcat('../../Results/ST/q_block_',num2str(i),'.pdf'),'-dpdf')


end
%%
% Loop through unblocked vessels
for i=[2, 18, 40]%1:length(unblock_ids)
figure(500+i); clf;
subplot(2,3,1); 
plot(t,p_c{unblock_ids(i),1}(:,1:2:end),'-r','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,2);
plot(t,p_hb{unblock_ids(i),1}(:,1:2:end),'-','Color',[0.82 0.52 0.63],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,3);
plot(t,p_h{unblock_ids(i),1}(:,1:2:end),'-','Color',[1.0 0.5 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,4); hold on;
plot(t,p_c{unblock_ids(i),2},'-b','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,5);
plot(t,p_hb{unblock_ids(i),2},'-','Color',[0.1 1.0 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,6);
plot(t,p_h{unblock_ids(i),2},'-','Color',[0.0 0.9 0.9],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

% print(strcat('../../Results/ST/p_unblock_',num2str(i),'.pdf'),'-dpdf')

%% Flow

figure(i+700); clf;
subplot(2,3,1); 
plot(t,q_c{unblock_ids(i),1}(:,1:2:end),'-r','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,2);
plot(t,q_hb{unblock_ids(i),1}(:,1:2:end),'-','Color',[0.82 0.52 0.63],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,3);
plot(t,q_h{unblock_ids(i),1}(:,1:2:end),'-','Color',[1.0 0.5 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,4); hold on;
plot(t,q_c{unblock_ids(i),2},'-b','LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,5);
plot(t,q_hb{unblock_ids(i),2},'-','Color',[0.1 1.0 0.8],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);

subplot(2,3,6);
plot(t,q_h{unblock_ids(i),2},'-','Color',[0.0 0.9 0.9],'LineWidth',2);
set(gca,'FontSize',20);
axis tight; %ylim([0 20]);
% print(strcat('../../Results/ST/q_unblock_',num2str(i),'.pdf'),'-dpdf')

end
%% Look at mean values throughout the alpha/beta branch
% Loop through blocked vessels
for i=[10,15, 25]%1:length(block_ids)%
figure(i+1000); clf;
subplot(2,2,1); hold on;
rac = sqrt(min(a_c{block_ids(i),1})./pi);
rbc = sqrt(min(a_c{block_ids(i),2})./pi);

rah = sqrt(min(a_h{block_ids(i),1})./pi);
rbh = sqrt(min(a_h{block_ids(i),2})./pi);

rahb = sqrt(min(a_hb{block_ids(i),1})./pi);
rbhb = sqrt(min(a_hb{block_ids(i),2})./pi);



plot(rac,mean(p_c{block_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,mean(p_hb{block_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,mean(p_h{block_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);
set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Alpha');

subplot(2,2,2); hold on;
plot(rbc,mean(p_c{block_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,mean(p_hb{block_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,mean(p_h{block_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Beta');

subplot(2,2,3); hold on;
plot(rac,mean(q_c{block_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,mean(q_hb{block_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,mean(q_h{block_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Flow - Alpha');

subplot(2,2,4); hold on;
plot(rbc,mean(q_c{block_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,mean(q_hb{block_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,mean(q_h{block_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);

title('Flow - Beta');


set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;

disp([i terminal(block_ids(i)) mean(q_c{block_ids(i),2}(:,1)) mean(q_hb{block_ids(i),2}(:,1))>mean(q_h{block_ids(i),2}(:,1))])
% print(strcat('../../Results/ST/meanblock_',num2str(i),'.eps'),'-depsc')
% pause(2)
end
%%
% Loop through unblocked vessels
for i=[2, 18, 40]%1:length(unblock_ids)%
figure(i+5000); clf;
subplot(2,2,1); hold on;
rac = sqrt(min(a_c{unblock_ids(i),1})./pi);
rbc = sqrt(min(a_c{unblock_ids(i),2})./pi);

rah = sqrt(min(a_h{unblock_ids(i),1})./pi);
rbh = sqrt(min(a_h{unblock_ids(i),2})./pi);

rahb = sqrt(min(a_hb{unblock_ids(i),1})./pi);
rbhb = sqrt(min(a_hb{unblock_ids(i),2})./pi);

plot(rac,mean(p_c{unblock_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,mean(p_hb{unblock_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,mean(p_h{unblock_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Alpha');

subplot(2,2,2); hold on;
plot(rbc,mean(p_c{unblock_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,mean(p_hb{unblock_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,mean(p_h{unblock_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);
set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Beta');

subplot(2,2,3); hold on;
plot(rac,mean(q_c{unblock_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,mean(q_hb{unblock_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,mean(q_h{unblock_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Flow - Alpha');

subplot(2,2,4); hold on;
plot(rbc,mean(q_c{unblock_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,mean(q_hb{unblock_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,mean(q_h{unblock_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);

title('Flow - Beta');


set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;

% print(strcat('../../Results/ST/meanunblock_',num2str(i),'.eps'),'-depsc')


pause (2)
end

%% Look at max pressure
% Loop through blocked vessels
for i=[10,15, 25]%1:length(block_ids)%
figure(i+10000); clf;
subplot(2,2,1); hold on;
rac = sqrt(min(a_c{block_ids(i),1})./pi);
rbc = sqrt(min(a_c{block_ids(i),2})./pi);

rah = sqrt(min(a_h{block_ids(i),1})./pi);
rbh = sqrt(min(a_h{block_ids(i),2})./pi);

rahb = sqrt(min(a_hb{block_ids(i),1})./pi);
rbhb = sqrt(min(a_hb{block_ids(i),2})./pi);



plot(rac,max(p_c{block_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,max(p_hb{block_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,max(p_h{block_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);
set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Alpha');

subplot(2,2,2); hold on;
plot(rbc,max(p_c{block_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,max(p_hb{block_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,max(p_h{block_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Beta');

subplot(2,2,3); hold on;
plot(rac,max(q_c{block_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,max(q_hb{block_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,max(q_h{block_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Flow - Alpha');

subplot(2,2,4); hold on;
plot(rbc,max(q_c{block_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,max(q_hb{block_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,max(q_h{block_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);

title('Flow - Beta');


set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;

disp([i terminal(block_ids(i)) mean(q_c{block_ids(i),2}(:,1)) mean(q_hb{block_ids(i),2}(:,1))>mean(q_h{block_ids(i),2}(:,1))])
% print(strcat('../../Results/ST/meanblock_',num2str(i),'.eps'),'-depsc')
% pause(2)
end
%%
% Loop through unblocked vessels
for i=[2, 18, 40]%1:length(unblock_ids)%
figure(i+50000); clf;
subplot(2,2,1); hold on;
rac = sqrt(min(a_c{unblock_ids(i),1})./pi);
rbc = sqrt(min(a_c{unblock_ids(i),2})./pi);

rah = sqrt(min(a_h{unblock_ids(i),1})./pi);
rbh = sqrt(min(a_h{unblock_ids(i),2})./pi);

rahb = sqrt(min(a_hb{unblock_ids(i),1})./pi);
rbhb = sqrt(min(a_hb{unblock_ids(i),2})./pi);

plot(rac,max(p_c{unblock_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,max(p_hb{unblock_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,max(p_h{unblock_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Alpha');

subplot(2,2,2); hold on;
plot(rbc,max(p_c{unblock_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,max(p_hb{unblock_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,max(p_h{unblock_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);
set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Pressure - Beta');

subplot(2,2,3); hold on;
plot(rac,max(q_c{unblock_ids(i),1}),'-o','MarkerSize',12,'Color','r','LineWidth',2);
plot(rahb,max(q_hb{unblock_ids(i),1}),'-.s','MarkerSize',12,'Color',[0.50 0 0],'LineWidth',2);
plot(rah,max(q_h{unblock_ids(i),1}),'-d','MarkerSize',12,'Color',[0.82 0.52 0.63],'LineWidth',2);

set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;
title('Flow - Alpha');

subplot(2,2,4); hold on;
plot(rbc,max(q_c{unblock_ids(i),2}),'-o','MarkerSize',12,'Color','b','LineWidth',2);
plot(rbhb,max(q_hb{unblock_ids(i),2}),'-.s','MarkerSize',12,'Color',[0.4824    0.4078    0.9333],'LineWidth',2);
plot(rbh,max(q_h{unblock_ids(i),2}),'-d','MarkerSize',12,'Color',[0 0.9 0.9 ],'LineWidth',2);

title('Flow - Beta');


set(gca,'FontSize',20,'XDir','reverse');
axis tight; grid on;

% print(strcat('../../Results/ST/meanunblock_',num2str(i),'.eps'),'-depsc')


pause (2)
end