close all; clear all; clc;

load CLdata.mat

figure;
gm = importGeometry('Healthy_Segmentation.stl');
pdegplot(gm,"FaceAlpha",0.3)

hold on

for i = 2:length(connectivity)
    vessel_info = cell2mat(vessel_details(i,2));
    radii = flip(vessel_info(:,4)./10); 
    distance = linspace(0, cell2mat(vessel_details(i,3))/10,length(radii));

    x = flip(vessel_details{i,2}(:,1));
    y = flip(vessel_details{i,2}(:,2));
    z = flip(vessel_details{i,2}(:,3));

    x_last = flip(vessel_details{i,2}(end,1));
    y_last = flip(vessel_details{i,2}(end,2));
    z_last = flip(vessel_details{i,2}(end,3));
    
    figure(1)
    hold on
    plot3(x,y,z,'linewidth',2,'color','k')
    plot3(x_last,y_last,z_last,'.','MarkerSize',30,'color','red')
    axis("off")
end


