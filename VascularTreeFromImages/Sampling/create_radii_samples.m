%%This runs the functions 'export_for_fluids_control' and
%%'export_for_fluids_sample' and creates the amount of samples you set the
%%index to

export_for_fluids_control('Healthy_menus.mat', 'connectivity_matrix.mat', 0.26)
for i = 1:10
    export_for_fluids_sample('Healthy_menus.mat', 'connectivity_matrix.mat', i, 0.26, 0.02)
end