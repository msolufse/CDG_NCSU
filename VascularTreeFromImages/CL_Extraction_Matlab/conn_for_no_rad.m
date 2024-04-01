function conn_for_no_rad(patientData)

    load(patientData);
    con_var = connectivity;
    
    % We want to initialize our connectivity matrix. We choose to initialize it
    % with zeros, as this will decrease the amount of work we have to do for
    % the terminal vessels.
    conn = zeros(height(con_var)-1, 4);
    conn_mat_str = table('Size',[height(con_var) 4], 'VariableTypes',{'cell','cell','cell','cell'});
    
    % The following lines creates a table version of our con_var matrix,
    % using the previous connectivity matrix, vessel_details, and vessel_radii.
    conn_mat_str(:, 1:4) = con_var(:, 1:4);
    conn_mat_str(:, 5) = vessel_details(:, 3);
    
    % This code converts the above table into a matrix of numerical values.
    % Note that the names of the vessels in conn_mat_str begin with a letter
    % and are followed by a number. We renamed the vessels by this number,
    % since there were no repeated numbers in the table of names above.
    conn(:, 1) = 0:height(con_var)-2;
    
    for n = 1:height(conn)
        conn(n, 2) = eraseBetween(string(table2array(conn_mat_str(n+1,2))), 1, 1);
        if string(table2array(conn_mat_str(n+1,2)))=="TERMINAL"
            conn(n,2) = 0;
            conn(n,3) = 0;
            conn(n,4) = 0;
            continue;
        end
    
        if isnan(cell2mat(table2array(conn_mat_str(n+1,3))))
            conn(n,3) = 0;
            conn(n,4) = 0;
            continue;
        else
            conn(n, 3) = double(eraseBetween(string(table2array(conn_mat_str(n+1,3))), 1, 1));
            
            if isempty(cell2mat(table2array(conn_mat_str(n+1,4))))
                conn(n,4) = 0;
            else
                conn(n, 4) = double(eraseBetween(string(table2array(conn_mat_str(n+1,4))), 1, 1));
            end
        end
    end
    
    term_ves = find(~conn(:,2)) - 1;
    
    save('connect.mat', "conn");
end
