%% This plots the pressure and flow for the MPA, LPA, and RPA

%%For pressure samples
%load time step
load pressure_mid_output/time_step.mat

%set index
for i = 1:10

    %Use append function to load matfiles within folder
    p = load(append("pressure_mid_output/pressure_mid",int2str(i),".mat"));
    figure(1)
    hold on
    
    %plot MPA pressure
    subplot(2,3,1); hold on
    plot(t,p.pressure_mid(:,1),'LineWidth',3,'Color','k')
    ylabel('Pressure mm/Hg')
    title('MPA')
    
    %plot LPA pressure
    subplot(2,3,2); hold on
    plot(t,p.pressure_mid(:,2),'LineWidth',3,'Color','k')
    title('LPA')
    
    %plot RPA pressure
    subplot(2,3,3); hold on
    title('RPA')
    plot(t,p.pressure_mid(:,128),'LineWidth',3,'Color','k')
end

%%For pressure control
%i should always be equal to 0
for i = 0

    %Use append function to load matfiles within folder
    p = load(append("pressure_mid_output/pressure_mid",int2str(i),".mat"));
    figure(1)
    hold on
    
    %plot MPA pressure
    subplot(2,3,1); hold on
    plot(t,p.pressure_mid(:,1),'LineWidth',3,'Color','r')
    
    %plot LPA pressure
    subplot(2,3,2); hold on
    plot(t,p.pressure_mid(:,2),'LineWidth',3,'Color','r')
    
    %plot RPA pressure
    subplot(2,3,3); hold on
    plot(t,p.pressure_mid(:,128),'LineWidth',3,'Color','r')
end

%%For flow samples
%load time step
load flow_mid_output/time_step.mat
for i = 1:10

    %Use append function to load matfiles within folder
    f = load(append("flow_mid_output/flow_mid",int2str(i),".mat"));
    figure(1)
    hold on
    
    %plot MPA pressure
    subplot(2,3,4); hold on
    plot(t,f.flow_mid(:,4),'LineWidth',3,'Color','k')
    ylabel('flow mm/Hg')
    xlabel('Time (s)')
    
    %plot LPA pressure
    subplot(2,3,5); hold on
    plot(t,f.flow_mid(:,5),'LineWidth',3,'Color','k')
    xlabel('Time (s)')

    %plot RPA pressure
    subplot(2,3,6); hold on
    plot(t,f.flow_mid(:,128),'LineWidth',3,'Color','k')
    xlabel('Time (s)')
end

%%For flow control
%i should always be equal to 0
for i = 0

%Use append function to load matfiles within folder
    f = load(append("flow_mid_output/flow_mid",int2str(i),".mat"));
    figure(1)
    hold on
    
    %plot MPA pressure
    subplot(2,3,4); hold on
    plot(t,f.flow_mid(:,4),'LineWidth',3,'Color','r')
    
    %plot LPA pressure
    subplot(2,3,5); hold on
    plot(t,f.flow_mid(:,5),'LineWidth',3,'Color','r')
    
    %plot RPA pressure
    subplot(2,3,6); hold on
    plot(t,f.flow_mid(:,128),'LineWidth',3,'Color','r')
end