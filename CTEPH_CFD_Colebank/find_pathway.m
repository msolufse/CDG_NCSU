% Take the terminal vessels and connectivity files and find the pathway
% from the MPA down to the terminal vessel. We just use this for
% preliminary results.

function paths = find_pathway(conn,terminal)

n = length(terminal);
paths = cell(n,1);
for i=1:n
    term_curr = terminal(i);
    path_curr = term_curr;
    id = find(conn(:,2)==term_curr);
    if isempty(id)
        id = find(conn(:,3)==term_curr);
    end
    while ~isempty(id)
        path_curr(end+1) = id;
        parent = conn(id,1);
        id = find(conn(:,2)==parent);
        if isempty(id)
            id = find(conn(:,3)==parent);
        end
    end
    paths{i} = path_curr;
end


end