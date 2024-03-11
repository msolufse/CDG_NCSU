%% This script will be used to understand the effects of geometric measurement on a
%  system of blood vessels, where we assume that the radius,length, and
%  connectivity are subject to uncertainty about their measurment and
%  population average. MJC 2/1/18
%%
function [tot_ves,tot_term,conn,dim_mat,terminal,geo3D] ...
    = create_data_fluids2(gen_ID,Stenosis,percent_sten,X_STEP)
% filename = 'Better_Lung_Data.mat';
% filename = 'Simvasc_arteries_Data.mat';
filename = 'new_pulm_Data.mat';
network = load(filename);
% X_STEP = 6;
connectivity = network.connectivity;
max_pts = 1024;
vessel_measurements = network.vessel_details;
% connectivity = network.connectivity;
%% Algorithm 7: Get simplified connectivity?
% [simple_connectivity,terminal] = find_simple_connectivity(connectivity, vessel_measurements);% Has number for parent, daughter1, and daughter2
[simple_connectivity,terminal] ...
    = find_simple_connectivity_V2(connectivity, vessel_measurements);
%% For the sake of connectivity check, we need to order the connectivity in a way that
% renumbers that connectivity in a better manner (and gives us bifurcations
% near each other. Try doing it by generation
% note: per generation, vessels increase (at most) by 2i. One vessel
% will hence go to 2*1 = 2 new vessels, then at generation 2 we get 4.
% Similarly, the size of the connectivity goes 1,3,7,15 (so up by 2 4 8
% ect).
generation = {};
generation{end+1} = simple_connectivity(1,1);
num_conn = 1;
generation_find = {1,2:3, 4:7, 8:15, 16:31, 32:63, 64:127, 128:255, 256:511, 512:1023, 1024:2047,2048:4095,4096:8191,8192:16383};
terminal_cell   = {};
terminal_cell{end+1} = 1;
new_daughters   = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,2048,4096,8192,16384];
for i=1:gen_ID+1
    if i~=1
        size_old = size(old_generation,1);
        num_conn = new_daughters(i-1); %Number of potential new daughters
        new_generation = zeros(size_old+num_conn,4);
        new_generation(1:size_old,:) = old_generation;
        new_inc = 1;
        clearID = [];
        for j=generation_find{i-1}
            d1 = old_generation(j,2);
            d2 = old_generation(j,3);
            d3 = old_generation(j,4);
            if d1~=0
                d1_id = find(simple_connectivity(:,1) == d1);
                d2_id = find(simple_connectivity(:,1) == d2);
                d3_id = find(simple_connectivity(:,1) == d3);
                if ~isempty(d1_id)
                    new_generation(size_old + 2*new_inc - 1,:) = simple_connectivity(d1_id,:);
                else
                    clearID(end+1) = size_old + 2*new_inc - 1;
                    %add terminal
                end
                if ~isempty(d2_id)
                    new_generation(size_old + 2*new_inc,:)  = simple_connectivity(d2_id,:);
                else
                    clearID(end+1) = size_old + 2*new_inc;
                    %add terminal
                end
                if ~isempty(d3_id)
                    new_generation(size_old + 2*new_inc,:)  = simple_connectivity(d3_id,:);
                else
                    clearID(end+1) = size_old + 2*new_inc;
                    %add terminal
                end
                new_inc = new_inc + 1;
            end
        end
        generation{end+1} = new_generation;
    else
        generation{end+1} = simple_connectivity(1,:);
    end
    old_generation = generation{end};
    new_terminal = [];
    for k=1:size(old_generation,1)
        check1 = old_generation(k,2);
        check2 = old_generation(k,3);
        check3 = old_generation(k,4);
        flag1  = isempty(find(old_generation(:,1) == check1,1));
        flag2  = isempty(find(old_generation(:,1) == check2,1));
        flag3  = isempty(find(old_generation(:,1) == check3,1));
        if flag1 == 1
            new_terminal(end+1) = check1;
        end
        if flag2 == 1
            new_terminal(end+1) = check2;
        end
        if flag3 == 1
            new_terminal(end+1) = check3;
        end
    end
    terminal_cell{end+1} = new_terminal;
end

for k=1:length(generation)
    gen_curr = generation{k};
    where = find(gen_curr(:,1) == 0);
    gen_curr(where,:) = [];
    generation{k} = gen_curr;
end

test_network = generation{gen_ID};
test_terminal = terminal_cell{gen_ID};
%% Define the connectivity
vessel_ids = unique(test_network,'stable')';
% Get rid of IDS that are zero
vessel_ids = vessel_ids(vessel_ids~=0);
num_vessels = max(size(unique(test_network)))-1; % Subtract one for zero indices
geo3d = cell(num_vessels+1,2);
% geo3d{1} = test_network;
%% Now define the connectivity for the model
measurements = zeros(num_vessels,2);
old_network = test_network;
old_terminal = test_terminal;
% Define the entry index for the vessels (i.e. vessel 300 -> entry 205)
abs_index = sort(vessel_ids);%[];
% for i=1:length(net3D)
%     abs_index(end+1) = net3D{i,2};
% end
% L = measurements(
for i=1:num_vessels %Go back through and reassign number (so you don't have 1 2 and 800 has your vessel #s)
    ves_curr = vessel_ids(i);
    ves_ID   = find(abs_index==ves_curr);
    [r_ves,c_ves] = find(old_network==ves_curr);
    for j = 1:length(r_ves)
        test_network(r_ves(j),c_ves(j)) = i;
    end
    [where_term] = find(old_terminal==ves_curr);
    for j = 1:length(where_term)
        test_terminal(where_term) = i;
    end
    %% CANNULA IS 0.86mm DIAMETER
    RL_precision = 3;
    measurements(i,1) = round(vessel_measurements{ves_curr,3}.*0.1,RL_precision);
    measurements(i,2) = round(vessel_measurements{ves_curr,4}.*0.1,RL_precision);
%     measurements(i,1) = round(dim_mat(ves_ID,2).*0.1,RL_precision);
%     measurements(i,2) = round(dim_mat(ves_ID,3).*0.1,RL_precision);
geo3D{i} = vessel_measurements{ves_curr,2};
end
test_terminal = test_terminal(test_terminal~=0);
test_terminal = sort(test_terminal);
%% THE MOST IMPORTANT STEP: Go back and fix indicies once the trifurcations have been dealt with
dlmwrite('connectivity.txt',test_network-1,'\t');
dlmwrite('terminal_vessels.txt',test_terminal-1,'\t');
%%
L       = zeros(num_vessels,1);
R       = zeros(num_vessels,1);
dim_mat = zeros(num_vessels,3);
id = 1;
for i=1:num_vessels
    L(i) = measurements(i,1);%round(vessel_measurements{vessels_21(i),3}.*0.1,3);
    R(i) = measurements(i,2);%round(vessel_measurements{vessels_21(i),4}.*0.1,3);
    %% VERY IMPORTANT TO INCLUDE THIS FOR SIMULATIONS
    % New here: Since we take X step sizes per cm, every length just
    % needs to be 1/X cm.
    if L(i) < (1/X_STEP)
        disp([i L(i) R(i)])
        L(i) = (2./X_STEP); %One point of spatial resolution
    end
    dim_mat(id,:) = round([L(i) R(i) R(i)],2);
    id = id+1;
end


tot_ves = max(max(test_network));
tot_term = length(test_terminal);

if test_network == 1
    conn = [1 0 0 0];
else
    conn = test_network;
end
terminal = test_terminal;
end