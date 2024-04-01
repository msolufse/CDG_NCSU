function standardize_networks(cm_name, num_vessels)

cm = load(cm_name);
conn_mat = cm.conn_mat;

new_conn = conn_mat;

% Go until the network has at most the desired number of vessels
while height(new_conn) > num_vessels

    radii        = sort(new_conn(:, 5));  
    min_radii    = radii(1);
    d_index      = find(new_conn(:, 5) == min_radii);
    not_terminal = new_conn(d_index, 2) ~= 0 || new_conn(d_index, 3) ~= 0 || new_conn(d_index, 4) ~= 0;

    % If the min radius vessel is not terminal, keep going until one is
    % found that is
    k = 2;
    if not_terminal
        while not_terminal
            % Iterate through radii vector
            min_radii    = radii(k);
            d_index      = find(new_conn(:, 5) == min_radii);
            not_terminal = new_conn(d_index, 2) ~= 0 || new_conn(d_index, 3) ~= 0 || new_conn(d_index, 4) ~= 0;
            k            = k+1;
        end
    end

    [p_index, ~] = find(new_conn(:, 2:4) == d_index-1); % Get parent index from found daughter

    d_1 = new_conn(p_index, 2);
    d_2 = new_conn(p_index, 3);
    d_3 = new_conn(p_index, 4);

    % Condition for if siblings are terminal
    if d_3 == 0 
        siblings_not_terminal = new_conn(d_1+1, 2) ~= 0 || new_conn(d_1+1, 3) ~= 0 || new_conn(d_1+1, 4) ~= 0 || new_conn(d_2+1, 2) ~= 0 || new_conn(d_2+1, 3) ~= 0 || new_conn(d_2+1, 4) ~= 0;
    else
        siblings_not_terminal = new_conn(d_1+1, 2) ~= 0 || new_conn(d_1+1, 3) ~= 0 || new_conn(d_1+1, 4) ~= 0 || new_conn(d_2+1, 2) ~= 0 || new_conn(d_2+1, 3) ~= 0 || new_conn(d_2+1, 4) ~= 0 || new_conn(d_3+1, 2) ~= 0 || new_conn(d_3+1, 3) ~= 0 || new_conn(d_3+1, 4) ~= 0;
    end

    % If siblings are not terminal vessels, we don't want to remove
    % generation, so keep going until a generation is found that is all
    % terminal
    if siblings_not_terminal
        while siblings_not_terminal

            min_radii    = radii(k);
            d_index      = find(new_conn(:, 5) == min_radii);
            not_terminal = new_conn(d_index, 2) ~= 0 || new_conn(d_index, 3) ~= 0 || new_conn(d_index, 4) ~= 0;
            k            = k+1;

            while not_terminal
                min_radii    = radii(k);
                d_index      = find(new_conn(:, 5) == min_radii);
                not_terminal = new_conn(d_index, 2) ~= 0 || new_conn(d_index, 3) ~= 0 || new_conn(d_index, 4) ~= 0;
                k            = k+1;
            end
        
            [p_index, ~] = find(new_conn(:, 2:4) == d_index-1);
        
            d_1 = new_conn(p_index, 2);
            d_2 = new_conn(p_index, 3);
            d_3 = new_conn(p_index, 4);
        
            if d_3 == 0
                siblings_not_terminal = new_conn(d_1+1, 2) ~= 0 || new_conn(d_1+1, 3) ~= 0 || new_conn(d_1+1, 4) ~= 0 || new_conn(d_2+1, 2) ~= 0 || new_conn(d_2+1, 3) ~= 0 || new_conn(d_2+1, 4) ~= 0;
            else
                siblings_not_terminal = new_conn(d_1+1, 2) ~= 0 || new_conn(d_1+1, 3) ~= 0 || new_conn(d_1+1, 4) ~= 0 || new_conn(d_2+1, 2) ~= 0 || new_conn(d_2+1, 3) ~= 0 || new_conn(d_2+1, 4) ~= 0 || new_conn(d_3+1, 2) ~= 0 || new_conn(d_3+1, 3) ~= 0 || new_conn(d_3+1, 4) ~= 0;
            end
            
        end
    end

    % Remove the daughters from network
    new_conn(d_1+1, :) = [];
    new_conn(find(new_conn(:, 1) == d_2), :) = [];
    if d_3 ~= 0
        new_conn(find(new_conn(:, 1) == d_3), :) = [];
    end

    new_conn(p_index, 2) = 0;
    new_conn(p_index, 3) = 0;
    new_conn(p_index, 4) = 0;

    % Update the connectivity
    for i = 1:height(new_conn)
        for j = 2:4
            old_vessel_index = new_conn(i, j);
            if old_vessel_index ~= 0
                new_vessel_index = find(new_conn(:, 1) == old_vessel_index);
                new_conn(i, j)   = new_vessel_index-1;
            end
        end
        new_conn(i, 1) = i-1;
    end

end

conn_mat = new_conn;

save(append(cm_name(1:end-4), '_standard.mat'), 'conn_mat')

end