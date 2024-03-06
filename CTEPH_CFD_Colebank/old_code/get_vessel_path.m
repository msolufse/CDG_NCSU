function out = get_vessel_path(conn,block,screen)
% Find the all the vessels distal to where blockages are located.
if ~isempty(block) && ~isempty(screen)
    block = [block screen];
elseif ~isempty(block)
    block = block;
elseif ~isempty(screen)
    block = screen;
else
    out = [];
    return;
end
nblock = length(block);
out = cell(nblock,1);
max_parent = max(conn(:,1));
for i=1:nblock
    path = [];
%     path(end+1) = block(i);
    d1 = conn(block(i),2);
    d2 = conn(block(i),3);
    d3 = conn(block(i),4);
    if d1>0
        path(end+1) = d1;
        path = recurr_ves_path(path,conn,d1,max_parent);
    end
    if d2>0
        path(end+1) = d2;
        path = recurr_ves_path(path,conn,d2,max_parent);
    end
    if d3>0
        path(end+1) = d3;
        path = recurr_ves_path(path,conn,d3,max_parent);
    end
    out{i} = path;
end

end


function path = recurr_ves_path(path,conn,i,max_parent)
if i>max_parent
    return;
end
d1 = conn(i,2);
d2 = conn(i,3);
d3 = conn(i,4);
if d1>0
    path(end+1) = d1;
    path = recurr_ves_path(path,conn,d1,max_parent);
end
if d2>0
    path(end+1) = d2;
    path = recurr_ves_path(path,conn,d2,max_parent);
end
if d3>0
    path(end+1) = d3;
    path = recurr_ves_path(path,conn,d3,max_parent);
end
end