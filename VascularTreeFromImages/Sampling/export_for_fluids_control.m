function export_for_fluids_control(vd_name, cm_name, min_length)

conn = prepare_for_fluids(vd_name, cm_name, min_length);

internal = [];
terminal = [];

for i = 1:height(conn)
    if conn(i, 2) == 0 && conn(i, 3) == 0 && conn(i, 4) == 0
        terminal = [terminal; conn(i, :)];
    else
        internal = [internal; conn(i, :)];
    end
end

for i = 1:height(internal)
    for j = 2:4
        index_in_internal = find(internal(:, 1) == internal(i, j));
        if isempty(index_in_internal)
            index_in_terminal = find(terminal(:, 1) == internal(i, j));
            internal(i, j) = index_in_terminal + height(internal) - 1;
        else
            internal(i, j) = index_in_internal - 1;
        end
    end
end

for i = 1:height(internal)
    internal(i, 1) = i-1;
end

for i = 1:height(terminal)
    terminal(i, 1) = i + height(internal) - 1;
end

parent_vessels = internal(:, 1:4);
terminal_vessels = terminal(:, 1)';

dimensions = [[internal(:, 6), internal(:, 5), internal(:, 5)]; [terminal(:, 6), terminal(:, 5), terminal(:, 5)]];
dimensions = ceil(dimensions*1000)/1000;

if ~isfolder('data_fluids/Sample0')
    mkdir('data_fluids/Sample0')
    addpath('data_fluids/Sample0')
end
    
dlmwrite('data_fluids/Sample0/connectivity.txt',parent_vessels, '\t');
dlmwrite('data_fluids/Sample0/terminal_vessels.txt',terminal_vessels, '\t');
dlmwrite('data_fluids/Sample0/dimensions.txt',dimensions, '\t');

end