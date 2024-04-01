function [simple_conn,terminal] = find_simple_connectivity_V2(conn_mat,details)
%% Note, this function processes based on the index of the matrix NOT the vessel number
M = size(conn_mat,1)-1;
% simple_conn = zeros(M,3);
temp_conn = zeros(M,4);
update_conn = zeros(M,4);
toss = [];
terminal = [];
for i=2:M+1
    clear who where radii_keep;
    k=2;
    radii_index = [];
    radii_comp  = [];
    daughters   = [];
    if any(toss == i) %Get rid of vessels that are part of trifurcation
        while ~isempty(conn_mat{i,k})% principal pathway
            if strcmp(conn_mat{i,2},'TERMINAL')
                break;
            end
            toss(end+1) = str2double(conn_mat{i,k}(2:end))+2;
            k=k+1;
        end
    elseif strcmp(conn_mat{i,2},'TERMINAL')
        terminal(end+1) = i;
        
    elseif isempty(conn_mat{i,4})
        temp_conn(i,1) = str2double(conn_mat{i,1}(2:end))+2;
        temp_conn(i,2) = str2double(conn_mat{i,2}(2:end))+2;
        temp_conn(i,3) = str2double(conn_mat{i,3}(2:end))+2;
    else
        if i==4 % This vessel is duplicated for some reason
            temp_conn(i,1) = str2double(conn_mat{i,1}(2:end))+2;
            temp_conn(i,2) = str2double(conn_mat{i,2}(2:end))+2;
            temp_conn(i,3) = str2double(conn_mat{i,4}(2:end))+2;
        else
            temp_conn(i,1) = str2double(conn_mat{i,1}(2:end))+2;
            temp_conn(i,2) = str2double(conn_mat{i,2}(2:end))+2;
            temp_conn(i,3) = str2double(conn_mat{i,3}(2:end))+2;
            temp_conn(i,4) = str2double(conn_mat{i,4}(2:end))+2;
        end
    end
end

temp_conn = temp_conn - update_conn;
parents = find(sum(temp_conn,2) > 0); %-1 for wanting vessel name (starts at 0)
simple_conn = temp_conn(parents,:); %+1 so indicies are greater than zero
end