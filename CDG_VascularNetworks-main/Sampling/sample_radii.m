function new_conn = sample_radii(conn,sample_number,min_radius)

new_conn = zeros(height(conn), width(conn));

new_conn(:, 1:4) = conn(:, 1:4);
new_conn(:, 6:7) = conn(:, 6:7);

for i = 1:height(conn)
    if i ~= 1
        [p_index, ~]  = find(conn(:, 2:4) == i-1);
        p_radius      = conn(p_index, 5);
        d_radius      = conn(i, 5);

        if conn(i, 7) ~= 0
            mu    = conn(i, 5);
            sigma = conn(i, 7);

            pd = makedist('normal', 'mu', mu, 'sigma', sigma);
        else
            scale = d_radius/p_radius;
            mu    = conn(p_index, 5)*scale;
            sigma = conn(p_index, 7)*scale;

            pd = makedist('normal', 'mu', mu, 'sigma', sigma);
        end

        radius_sample = random(pd);
        % parent_sample = new_conn(p_index, 5);

        % while radius_sample >= parent_sample
        %     radius_sample = random(pd);
        % end

        while radius_sample <= min_radius %.021
            radius_sample = random(pd);
        end

    else
        mu    = conn(i, 5);
        sigma = conn(i, 7);

        pd = makedist('normal', 'mu', mu, 'sigma', sigma);

        radius_sample = random(pd);
    end

    new_conn(i, 5) = radius_sample;
end

conn_mat = new_conn;

save(append('connectivity_matrices/conn_mat_sample', num2str(sample_number), '.mat'), 'conn_mat')

end