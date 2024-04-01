function conn_mat = prepare_for_fluids(vd_name, cm_name, min_length)

vd = load(vd_name);
cm = load(cm_name);

vessel_details = vd.vessel_details;
connectivity = vd.connectivity;

conn_mat = cm.conn_mat;

vessel_length = 0;
for i = 1:height(conn_mat)
    vessel_length = conn_mat(i, 6);
    if vessel_length <= min_length
        vessel_length = min_length;
        vessel_details{i+1, 3} = vessel_length;
        conn_mat(i, 6) = vessel_length;
    end
end

for i = 2:height(vessel_details)
    vessel_details{i, 3} = vessel_details{i, 3}*10;
    vessel_details{i, 4} = vessel_details{i, 4}*10;
end

save(append(vd_name(1:end-4), '_updated.mat'), 'vessel_details', 'connectivity');
save(append(cm_name(1:end-4), '_updated.mat'), 'conn_mat');

end