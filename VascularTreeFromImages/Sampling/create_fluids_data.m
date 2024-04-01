function create_fluids_data(conn_mat, sample_number)

conn = conn_mat;

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

if ~isfolder(append('data_fluids/Sample', num2str(sample_number)))
    mkdir(append('data_fluids/Sample', num2str(sample_number)))
    addpath(append('data_fluids/Sample', num2str(sample_number)))
end

dlmwrite(append('data_fluids/Sample', num2str(sample_number), '/connectivity.txt'),parent_vessels, '\t');
dlmwrite(append('data_fluids/Sample', num2str(sample_number), '/terminal_vessels.txt'),terminal_vessels, '\t');
dlmwrite(append('data_fluids/Sample', num2str(sample_number), '/dimensions.txt'),dimensions, '\t');

end