% Find which terminal vessels are affected by stenoses to make plots within
% the structured tree
clear; clc; close all;
%%
load ../term_conn_sten.mat
% Get path where stenosis affects
path = get_vessel_path(conn,block,screen);


% Decide which vessels to use
block_ids   = [];
unblock_ids = 1:length(terminal);
for path_id = 1:size(path,1)
    for i=1:length(path{path_id})
        if any(path{path_id}(i)==terminal)
            block_ids(end+1) = find(path{path_id}(i)==terminal);
        end
    end
end
block_ids = unique(block_ids);
unblock_ids(block_ids) = [];
save('term_IDs','block_ids','unblock_ids');
