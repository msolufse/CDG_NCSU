%% Recusive Function for finding networks
function [vessels_Org, new_network, conn_mat, name, ID, connID, charID] =  find_network_cl(parent, ...
    old_network, vessels_Org, conn_mat, name, ID, connID, charID, startnodes)
new_network = old_network;
%% INITIAL COST SHOULD BE THE DISTANCE FROM THE END OF THE DAUGHTER VESSEL TO ITS PARENT
% This is an area of the code that can often fail, since we don't account
% for how vessels might be curved
%% This is an attempt to change the way we scan for where the bifurcation points should line up with centerlines
NetworkSize = max(size(old_network));
%% As a manual check for where the bifurcations are showing up
% figure(200);
if name-1 >= size(vessels_Org,1)
    display(name);
end
if strcmp(char(charID),'B')
    1
end
parentnode = parent(1,1:3);

% Test to see if its terminal
terminalflag = 1;
for i=1:length(startnodes)
   if parentnode == startnodes(i,:)
       terminalflag = 0;
       break
   end
end
if terminalflag == 1
       conn_mat{connID,2} = 'TERMINAL';
    return
end

%% Find the bifurcation vector of interest
where = 0;
daughtercounter = 0;
daughterindex = [];
for i=1:NetworkSize
    if ~isempty(new_network{i}) && sum(parentnode==new_network{i}(end,1:3))==3
            daughtercounter = daughtercounter + 1;
            daughterindex(end+1) = i;
    end
end
%% Added loop to address any problems
if daughtercounter == 1
    warning('Only one daughter branch found. Attempting to other daughter.\n')
end
if daughtercounter == 0
       conn_mat{connID,2} = 'TERMINAL';
    return
end
ID_inner = ID+1;
conn_inner = connID;
charIDinner = charID + 1;
% As a note, here are the following uses for each increment:
% connID: denotes the row we are appending in connectivity matrix
% conn_inner: denotes ".." in RECURSIVE connectivity matrix
% name:  describes the name of the vessel (i.e. V1, V2, ect).
% bif_inc: keeps track of how many daughters we have recursively assesssed
% ID: keeps track of which vessel in the new_network we are on
bif_inc = 1;
while bif_inc <= daughtercounter % Figure out how many new vessels we attach
    if daughtercounter==3
        2
    end
    name = name+1; % Update which vessel we are at
    ID = ID + 1;% New spot in organized vessel array
    conn_inner = conn_inner+1; %Next row after parent
    conn_mat{connID, bif_inc+1} = strcat(char(charID),num2str(name));%Connectivity of parent
    conn_mat{connID,9} = [conn_mat{connID,9} name];
    vessels_Org{ID,1} = strcat(char(charID),num2str(name));
    vessels_Org{ID,2} = old_network{daughterindex(bif_inc)};% Based on the lowest cost
    new_network{daughterindex(bif_inc)} = [];  % values for bifurcation - vessel
    startD = vessels_Org{ID,2}(:,1:3);
    conn_mat{conn_inner,1} = vessels_Org{ID,1};
    %%
    [vessels_Org, new_network, conn_mat, name, ID, conn_inner, charIDinner] =  find_network_cl(startD, ...
        new_network, vessels_Org, conn_mat, name, ID, conn_inner, charIDinner, startnodes);
    bif_inc = bif_inc + 1;
end
connID = conn_inner;
end
