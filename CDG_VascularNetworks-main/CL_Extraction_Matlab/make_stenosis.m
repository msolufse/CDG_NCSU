%%
% Function that will take a geometry in and change everything to have a
% stenosis in one of the vessels.

function [conn, term, dim_mat, sten_loc] = make_stenosis(conn,term,dim_mat,...
                        sten_loc)

% Now we want to update the connectivity based on the stenosis location
% conn = conn;% + (conn>sten_loc);
% term = term;% + (term>sten_loc);

% Loop through all the screens

if ~isempty(sten_loc) % If there are stenoses, account for them
    [ids,~,~] = find(conn(:,1)==sten_loc);
    conn(ids,4) = -1;
end



end