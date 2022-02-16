%% This script will be used to understand the effects of geometric measurement on a
%  system of blood vessels, where we assume that the radius,length, and
%  connectivity are subject to uncertainty about their measurment and
%  population average. MJC 2/1/18
%% 'ves_conn','term','dim_mat'
function [geo3d,tot_ves,tot_term] = create_data_fluids(gen_ID, ves_conn, term, dim_mat,net3D,BIFSCALE)%,RSTD)
new_connectivity = ves_conn;
%% For the sake of connectivity check, we need to order the connectivity in a way that
% renumbers that connectivity in a better manner (and gives us bifurcations
% near each other. Try doing it by generation
% note: per generation, vessels increase (at most) by 2i. One vessel
% will hence go to 2*1 = 2 new vessels, then at generation 2 we get 4.
% Similarly, the size of the connectivity goes 1,3,7,15 (so up by 2 4 8
% ect).
generation = {};
generation{end+1} = ves_conn(1,1);
num_conn = 1;
generation_find = {1,2:3, 4:7, 8:15, 16:31, 32:63, 64:127, 128:255, 256:511, 512:1023, 1024:2047,2048:4095,4096:8191};
terminal_cell   = {};
terminal_cell{end+1} = 1;
new_daughters   = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,2048,4096,8192];
for i=1:gen_ID+1
    if i~=1
        size_old = size(old_generation,1);
        num_conn = new_daughters(i-1); %Number of potential new daughters
        new_generation = zeros(size_old+num_conn,3);
        new_generation(1:size_old,:) = old_generation;
        new_inc = 1;
        clearID = [];
        for j=generation_find{i-1}
            d1 = old_generation(j,2);
            d2 = old_generation(j,3);
            if d1~=0
                d1_id = find(new_connectivity(:,1) == d1);
                d2_id = find(new_connectivity(:,1) == d2);
                if ~isempty(d1_id)
                    new_generation(size_old + 2*new_inc - 1,:) = new_connectivity(d1_id,:);
                else
                    clearID(end+1) = size_old + 2*new_inc - 1;
                    %add terminal
                end
                if ~isempty(d2_id)
                    new_generation(size_old + 2*new_inc,:)  = new_connectivity(d2_id,:);
                else
                    clearID(end+1) = size_old + 2*new_inc;
                    %add terminal
                end
                new_inc = new_inc + 1;
            end
        end
        generation{end+1} = new_generation;
    else
        generation{end+1} = new_connectivity(1,:);
    end
    old_generation = generation{end};
    new_terminal = [];
    for k=1:size(old_generation,1)
        check1 = old_generation(k,2);
        check2 = old_generation(k,3);
        flag1  = isempty(find(old_generation(:,1) == check1,1));
        flag2  = isempty(find(old_generation(:,1) == check2,1));
        if flag1 == 1
            new_terminal(end+1) = check1;
        end
        if flag2 == 1
            new_terminal(end+1) = check2;
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
num_vessels = max(size(unique(test_network)));%size(vessel_measurements,1)-1;%length(vessels_21);
geo3d = cell(num_vessels+1,2);
% geo3d{1} = test_network;
%% Now define the connectivity for the model
measurements = zeros(num_vessels,2);
old_network = test_network;
old_terminal = test_terminal;
% Define the entry index for the vessels (i.e. vessel 300 -> entry 205)
abs_index = [];
for i=1:length(net3D)
    abs_index(end+1) = net3D{i,2};
end
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
    if i==1
        % FOR HARD CODING MPA
        if ~isempty(BIFSCALE)
            dim_mat(1,3) = BIFSCALE(1);
            dim_mat(2,3) = BIFSCALE(2);
            dim_mat(3,3) = BIFSCALE(3);
        end
        alpha = 0.43./dim_mat(ves_ID,3);
        dim_mat = dim_mat.*alpha;
%         RSTD    = RSTD.*alpha*0.1;
    end
    measurements(i,1) = round(dim_mat(ves_ID,2).*0.1,RL_precision);
    measurements(i,2) = round(dim_mat(ves_ID,3).*0.1,RL_precision);
end

%% Added by MJC 3/11/2020
% Before saving this geometry, we need to verify that all the vessels
% satisfy the area condition, Ap>Ad1,Ad2
% p  = test_network(:,1);
% d1 = test_network(:,2);
% d2 = test_network(:,3);
% 
% figure(989); clf; 
% hold on; 
% plot(measurements(p,2),'k*','MarkerSize',10); 
% plot(measurements(d1,2),'ro','MarkerSize',10); 
% plot(measurements(d2,2),'co','MarkerSize',10);
% for i=length(p):-1:1
%    r_i = measurements([p(i) d1(i) d2(i)],2);
%    if r_i(1) < r_i(2)
%        measurements(p(i),2) = r_i(2);
%    end
%    if r_i(1) < r_i(3)
%        measurements(p(i),2) = r_i(3);
%    end
% end
% plot(measurements(p,2),'ks','MarkerSize',10); 
% plot(measurements(d1,2),'r*','MarkerSize',10); 
% plot(measurements(d2,2),'c*','MarkerSize',10);
% 2


%% PLOT THE 3D NETWORK USED

figure; hold on;
for i=1:num_vessels
    plot_id = find(vessel_ids(i)==abs_index);
    plot3(net3D{plot_id,1}(:,1),net3D{plot_id,1}(:,2),net3D{plot_id,1}(:,3),'o');
    text(net3D{plot_id,1}(floor(end/2),1),net3D{plot_id,1}(floor(end/2),2),net3D{plot_id,1}(floor(end/2),3),strcat('V',num2str(i)),'FontSize',20);
           geo3d{i,1} = net3D{vessel_ids(i)};
           geo3d{i,2} = measurements(i,2);
end

%% THE MOST IMPORTANT STEP: Go back and fix indicies once the trifurcations have been dealt with
dlmwrite('connectivity.txt',test_network-1,'\t');
dlmwrite('terminal_vessels.txt',test_terminal-1,'\t');
%%
dim_mat2 = zeros(num_vessels,2);
id = 1;
for i=1:num_vessels
    L(i) = measurements(i,1);%round(vessel_measurements{vessels_21(i),3}.*0.1,3);
    R(i) = measurements(i,2);%round(vessel_measurements{vessels_21(i),4}.*0.1,3);
    %% VERY IMPORTANT TO INCLUDE THIS FOR SIMULATIONS
    % Except we are sampling, so this doesn't really need to be applied
    x_spat = 1./40;
    if L(i) < x_spat
%         disp([i L(i) R(i)])
%         L(i) = 1./20;
        L(i) = x_spat;%4.*R(i);
    end
    dim_mat2(id,:) = [L(i) R(i)];
    id = id+1;
end

temp_network = test_network; %Make a temporary network that will be reduced
temp_terminal = test_terminal;%Make a temporary terminal array that will be reduced

%%
dim_mat2 = round([L' R'],RL_precision);
%%
dlmwrite('Dimensions.txt',dim_mat2,'\t')
%% Add a loop here to do the following: we want to reorder the
% connectivity matrix and the terminal vessels so we get the
% largest values first (i.e. vessel 27 before 1) so the C code can
% parse it correctly.
temp_terminal = sort(temp_terminal);
[who,where] = sort(temp_network(:,1));
temp_network = temp_network(where,:);
dlmwrite('connectivity.txt',temp_network-1,'\t');
dlmwrite('terminal_vessels.txt',temp_terminal-1,'\t');

tot_ves = max(max(temp_network));
tot_term = length(temp_terminal);

end