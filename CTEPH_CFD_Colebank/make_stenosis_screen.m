%%
% Function that will take a geometry in and change everything to have a
% stenosis in one of the vessels.

function [conn, term, dim_mat, sten_loc,screen_loc] = make_stenosis_screen(conn,term,dim_mat,...
                        sten_loc,screen_loc,screen_L,num_pts)

% Now we want to update the connectivity based on the stenosis location
% conn = conn;% + (conn>sten_loc);
% term = term;% + (term>sten_loc);

% Loop through all the screens
n_sten = length(screen_loc);

if ~isempty(sten_loc) % If there are stenoses, account for them
    [ids,~,~] = find(conn(:,1)==sten_loc);
    conn(ids,4) = -1;
end

for i=1:n_sten
    max_ves = max(conn(:));
    max_parent = max(conn(:,1));

    if screen_loc(i)<=max_parent
        % Find where the parent is
        p_id = find(conn(:,1)==screen_loc(i));
        % find the minimum of the daughters (for indexing purposes)
        % If we are at a trifurcation, consider the 4th vessel
        if conn(p_id,4)>0
            min_d = min(conn(p_id,2),min(conn(p_id,3),conn(p_id,4)));
        else
            min_d = min(conn(p_id,2),conn(p_id,3));
        end
        %redefine the connectivity
        ids = conn>screen_loc(i); % Find all the vessels greater than the screen
        conn = conn+ids;          % Increment them by 1
        temp = conn(1:p_id,:);    % Define the temporary connectivity that will be updated
        temp(end,:) = [p_id p_id+1 0 0]; % Set up screen vessel
        temp(end+1,:) = [conn(p_id,1)+1 conn(p_id,2) conn(p_id,3) conn(p_id,4)];
        temp(end+1:end+(max_parent-p_id),:) = conn(p_id+1:end,:);
        conn = temp;
        ids = term>screen_loc(i);
        term = term+ids;
        
        if any(sten_loc>screen_loc(i))
            bool = sten_loc>screen_loc(i);
            sten_loc = sten_loc+bool;
        end
        bool = screen_loc>screen_loc(i);
        screen_loc = screen_loc+bool;
        
        % Update the radius and length
        L = dim_mat(screen_loc(i),1);
        R = dim_mat(screen_loc(i),2);
        
        L_ob = L./2 - (screen_L.*L)./2; % Stenosis at halfway
        L_ob = max(L_ob,1/num_pts);    % Stenosis at halfway
        
%         L_ob = L - (sten_L.*L)./2;  % Stenosis at bif
%         L_ob = max(L_ob,1/num_pts);
        
        temp = dim_mat(1:screen_loc(i),:);
        temp(end,1) = L_ob;
        temp(end+1,:) = temp(end,:);
        temp(end+1:end+(max_ves-p_id),:) = dim_mat(p_id+1:end,:);
        dim_mat = temp;
%         screen_loc = screen_loc+1;
    else
        conn(end+1,:) = [screen_loc(i) max_ves+1 0 0];
        term(term==screen_loc(i)) = max_ves+1;
        term = sort(term); % Need to resort so terminal vessels are in decreasing order.
        
        % Update the radius & length
        L = dim_mat(screen_loc(i),1);
        R = dim_mat(screen_loc(i),2);
        
        L_ob = L./2 - (screen_L.*L)./2;
        L_ob = max(L_ob,1/num_pts);
        dim_mat(screen_loc(i),1) = L_ob;
        dim_mat(end+1,:) = [L_ob,R,R];
    end
    
end



end