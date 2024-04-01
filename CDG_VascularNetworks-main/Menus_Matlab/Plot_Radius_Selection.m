load Patients/Healthy_Patient.mat;
load Data/vessel_radii.mat;
load Data/changepoint_locations.mat;


for i = 3%:length(vessel_radii)

    l = linspace(0, vessel_details{i, 3}, height(vessel_details{i, 2}))';
    r = vessel_details{i, 2}(1:end, 4);

    radius_estimation = vessel_radii{i, 2}.*ones(height(vessel_details{i, 2}), 1);
    radius_plus = (vessel_radii{i, 2} + vessel_radii{i, 3}).*ones(height(vessel_details{i, 2}), 1);
    radius_minus = (vessel_radii{i, 2} - vessel_radii{i, 3}).*ones(height(vessel_details{i, 2}), 1);

    f = figure(i-1);
    hold on
    p1 = plot(l, (r), '.', 'MarkerSize', 20, 'LineWidth', 2);
    p2 = plot(vessel_radii{i, 4}(:, 1), vessel_radii{i, 4}(:, 2), '.', 'color', 'green', 'MarkerSize', 20, 'LineWidth', 2);

    p3 = plot(l, radius_estimation, 'color', 'black', 'LineWidth', 2.0);
    p4 = plot(l, radius_plus, '--', 'color', 'black', 'LineWidth', 2.0);
    plot(l, radius_minus, '--', 'color', 'black', 'LineWidth', 2.0)
    cp_l = changepoint_locations{i,2}(:,2);
    cp_r = changepoint_locations{i,2}(:,1);
    plot(cp_l, cp_r, '.', 'MarkerSize', 20,'Color','red')

    grid minor
    xlabel('Length (cm)')
    ylabel('Radius (cm)')
    set(gca, 'fontsize', 20)
    
    if vessel_details{i, 11} == "TERMINAL,,,,,"
        title("Vessel " + vessel_details{i, 1} + " (Terminal)")
    else
        title("Vessel " + vessel_details{i, 1})
    end
    hold off

    legend([p1 p2 p3 p4], {'Observed Radii', 'Region Considered', 'Est. Radius', ...
        'Error'}, 'Location', 'best', 'FontSize', 12)
    print(['Figures/Radius_Region' int2str(i-1)],'-dpng')
end