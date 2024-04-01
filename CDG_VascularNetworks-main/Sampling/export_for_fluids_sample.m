function export_for_fluids_sample(vd_name, cm_name, sample_num, min_length, min_radius)

conn_mat = prepare_for_fluids(vd_name, cm_name, min_length);
sample_conn = sample_radii(conn_mat, sample_num, min_radius);
create_fluids_data(sample_conn, sample_num)

end