%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find_nodes.m File
% Author: M. J. Colebank
% Last update: 4/20/19
%
% Description: This code takes in the VMTK files (which contain repeated
% centerline points) and separates the file into unique points. Once these
% unique points are found, the nodes/arcs (corresponding to bifurcations
% and vessels) are determined, and can then be passed to the network
% connectivity finder.
%
% nodes: the junctions points at a given bifurcations
%
% split_arcs: the arcs for individual, unique vessels
%
% last_point: the lowest point, which allows for adjusting the tree
% structure (if that is something required)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes, split_arcs, last_point] = find_nodes(cl_unswapped)
cl = cl_unswapped(:,[2 3 4 1]);
[uni_mat, ~, unique_location] = unique(cl,'rows','stable');
N2 = size(uni_mat,1);
count = zeros(N2,1);
for i=1:N2 %only loop through unique points
    where = [];
    findings = (i==unique_location);
    wherefindings = find(i == unique_location);
    diff_find = diff(wherefindings);
    if any(diff_find == 1) && length(diff_find) > 1
        where = find(diff_find == 1); %Delete any duplicates in the array
        cl(wherefindings(where+1),:) = [];
        [~, ~, unique_location] = unique(cl,'rows','stable');
    end
    count(i) = sum(findings) - length(where);
end
[uni_mat, ~, unique_location] = unique(cl,'rows','stable');
% lowest_point_id = find(uni_mat(:,3) ==min(uni_mat(:,3)));
lowest_point = cl(end,:); % FOR NOW, JUST TAKE THE LAST POINT IN THE CL FILE (MJC 2/25/2020)
lowest_vec   = cl(end,:)-cl(end-1,:);
%% Find Segment ID's (jump discontinuities)
N = size(cl,1);
N2 = size(uni_mat,1);
r = cl(:,4);
xyz = cl(:,1:3);
vectors = zeros(N,5);
M = N-1;
for i=1:M
    vectors(i,1:3) = [xyz(i+1,1) - xyz(i,1),...
        xyz(i+1,2) - xyz(i,2), xyz(i+1,3) - xyz(i,3)];
    vectors(i,4) = sqrt( (xyz(i+1,1) - xyz(i,1)).^2 ...
        + (xyz(i+1,2) - xyz(i,2)).^2 + (xyz(i+1,3) - xyz(i,3)).^2);
    vectors(i,5) = r(i+1)+r(i);
end

% Now compare the distance between two (x,y,z) points and their total
% radius distance between. If its greater than 2r, its a jump in the CL.
dist_r = vectors(:,4:5);
segmentID = [1 (find(diff(dist_r')<0)) + 1 N];


% Here is a plotting routine to see if you've found all the
% new starts to full pathways.
figure(1); clf; hold on;
plot3(uni_mat(:,1),uni_mat(:,2),uni_mat(:,3),'o');
plot3(xyz(segmentID,1),xyz(segmentID,2), ...
    xyz(segmentID,3),'ks','MarkerSize',10,'MarkerFace','k');
%% ADDED TO ADDRESS THE PROBLEM OF NOT FINDING THE ROOT VESSEL
prompt = 'Was the root vessel found? 0 - No, 1- Yes\n';
flag = input(prompt);
% If the root vessel was not found (as can happen for short root arteries),
% then we simply take the root vessel as an extension of the lowest z-axis
% point in the tree.
if flag == 0
    %    new_root_z = lowest_point(3);
    %    xroot      = lowest_point(1)*ones(10,1);
    %    yroot      = lowest_point(2)*ones(10,1);
    %    zroot      = linspace(0.99*new_root_z,new_root_z.*(1.*lowest_point(4)),10)';
    %    rroot      = lowest_point(4)*ones(10,1);
    %    new_root   = [xroot yroot zroot rroot];
    %    lowest_point = new_root(1,:);
    %    cl(end+1:end+10,:) = new_root;
    % Use the difference vector to determine which to put the new root of
    % the vessel.
    xroot_end = 10*lowest_vec(1)+lowest_point(1);
    yroot_end = 10*lowest_vec(2)+lowest_point(2);
    zroot_end = 10*lowest_vec(3)+lowest_point(3);
    
    xroot = linspace(lowest_point(1),xroot_end,10);
    yroot = linspace(lowest_point(2),yroot_end,10);
    zroot = linspace(lowest_point(3),zroot_end,10);
    rroot      = lowest_point(4)*ones(10,1);
    new_root   = [xroot' yroot' zroot' rroot];
    cl(end+1:end+10,:) = new_root;
end
%% Need a loop here to see if we have endpoints appended to each vessel.
% If this is the case, we need to append the last points to each vessel.
if mode(diff(segmentID)) == 1
    endofnewvessels = find(diff(segmentID) == 1);
    endofnewvessels(end+1) = endofnewvessels(end)+1; %Append an N because we previously appended N to the end (i.e. the difference came out to be zero, so we didn't detect it)
    % We go from 2:end here because the first entry in endofnewvessels is
    % actually the end of a full segment
    A = cl(segmentID(endofnewvessels(1:end)),:);
    cl_new = [];
    for i=1:size(A,1)./2 %A should always be of even dimensions
        length_interval = length(segmentID(i):segmentID(i+1)-1);
        cl_new(end+1,:) = A((2*i)-1,:);
        cl_new(end+1:end+length_interval,:) =  cl(segmentID(i):segmentID(i+1)-1,:);
    end
    %% Now recalculate the segmentID jumps
    cl = cl_new;
    [uni_mat, ~, unique_location] = unique(cl,'rows','stable');
    lowest_point_id = find(uni_mat(:,3) ==min(uni_mat(:,3)));
    lowest_point = uni_mat(lowest_point_id,:);
    N = size(cl,1);
    N2 = size(uni_mat,1);
    r = cl(:,4);
    xyz = cl(:,1:3);
    vectors = zeros(N,5);
    M = N-1;
    for i=1:M
        vectors(i,1:3) = [xyz(i+1,1) - xyz(i,1),...
            xyz(i+1,2) - xyz(i,2), xyz(i+1,3) - xyz(i,3)];
        vectors(i,4) = sqrt( (xyz(i+1,1) - xyz(i,1)).^2 ...
            + (xyz(i+1,2) - xyz(i,2)).^2 + (xyz(i+1,3) - xyz(i,3)).^2);
        vectors(i,5) = r(i+1)+r(i);
    end
    % Now compare the distance between two (x,y,z) points and their total
    % radius distance between. If its greater than 2r, its a jump in the CL.
    dist_r = vectors(:,4:5);
    segmentID = [1 (find(diff(dist_r')<0)) + 1 N];
end
% Now we want to take the starts of each centerline path that is located in
% the full CL file, and map it to a point in the unique matrix.
uni_segmentID = unique_location(segmentID);
uni_segmentID(end) = size(uni_mat,1);

%% Try and sort through the entire centerline file and break full CL paths
% into individual segments. Save the nodes inbetween and compare later.
% Note that this should ONLY USE THE CL FILE AND SEGMENT_ID.
len_segID = 1:(length(segmentID)-1);
temp_nodes = []; %might not need this in the long run
final_nodes = [];
split_arcs = {};
start_nodes = []; end_nodes = [];
last_point = cl(N,:);
ves_node_ID = []; %%Keeps track of the vesID,and node IDs
for i = len_segID
    final_nodes(end+1,:) = cl(segmentID(i),:); %Append the start of each full centerline path
    cl_curr = cl(segmentID(i):segmentID(i+1)-1,:); %Centerline to compare
    start_nodes(end+1,:) = cl_curr(1,:);
    not_seg = len_segID(len_segID~=i);
    len_curr = size(cl_curr,1);
    temp_intersect = []; %To save the place at which there is an intersection
    for j=not_seg
        cl_comp = cl(segmentID(j):segmentID(j+1)-1,:);
        len_comp = size(cl_comp,1);
        ones_mat = ones(len_comp,4);
        for k=1:len_curr
            comp = sum(abs(ones_mat.*cl_curr(k,:) - cl_comp),2);
            where = find(comp == 0); %Find identical points
            if ~isempty(where)
                temp_intersect(end+1) = where(1)+ segmentID(j); %Add segmentID since we are looking at a small segment of CL file
                break
            end
        end
    end
    intersect_nodes = unique(cl(temp_intersect,:),'rows');
    intersect_nodes(end+1,:) = last_point;
    len_intersect = size(intersect_nodes,1);
    temp_nodes(end+1:end+len_intersect,:) = intersect_nodes;
    %% Before splitting the arcs, we need to findout what nodes have already been used.
    %  Then go through the original arc and separate the subarcs
    start_subarc = 1;
    temp_start = start_nodes(end,:); %Always assign to the last value, since we appended it at the beginning of the loop
    for t=1:len_curr
        if any(cl_curr(t,:) == intersect_nodes)
            if ~isempty(end_nodes)
                temp_end = cl_curr(t,:);
                dist_start = abs(sum(ones(size(start_nodes,1),4).*temp_start - start_nodes,2));
                dist_end   = abs(sum(ones(size(end_nodes,1),4).*temp_end - end_nodes,2));
                %%
                % Below, we don't include the last point in dist_start,
                % since the previous vessel's starting point is located at
                % the end of this array.
                
                if size(dist_start,1) == 1
                    end_nodes(end+1,:) = cl_curr(t,:);
                    start_nodes(end+1,:) = cl_curr(t,:);
                    split_arcs{end+1} = cl_curr(start_subarc:t,:);
                    ves_node_ID(end+1,1:3) = [size(split_arcs,2) size(start_nodes,1) size(end_nodes,1)];
                    start_subarc = t;
                    temp_start = cl_curr(t,:);
                elseif min(dist_start(1:end-1)) == 0 && min(dist_end) == 0 %If both nodes have been found before; don't include this subarc
                elseif min(dist_start(1:end-1))~=0 && min(dist_end) == 0 % Only the start point is new
                    split_arcs{end+1} = cl_curr(start_subarc:t,:);
                    min_start = find(min(dist_start(1:end-1)) == dist_start);
                    min_end   = find(min(dist_end(1:end-1)) == dist_end);
                    ves_node_ID(end+1,1:3) = [size(split_arcs,2) min_start min_end];
                    start_subarc = t;
                    temp_start = cl_curr(t,:);
                else
                    end_nodes(end+1,:) = cl_curr(t,:);
                    if sum(cl_curr(t,:) ~= lowest_point)==4
                        start_nodes(end+1,:) = cl_curr(t,:);
                    end
                    split_arcs{end+1} = cl_curr(start_subarc:t,:);
                    ves_node_ID(end+1,1:3) = [size(split_arcs,2) size(start_nodes,1) size(end_nodes,1)];
                    start_subarc = t;
                    temp_start = cl_curr(t,:);
                end
            else
                end_nodes(end+1,:) = cl_curr(t,:);
                split_arcs{end+1} = cl_curr(start_subarc:t,:);
                ves_node_ID(end+1,1:3) = [size(split_arcs,2) size(start_nodes,1) size(end_nodes,1)];
                start_subarc = t;
                temp_start = cl_curr(t,:);
            end
            %                        plot3(split_arcs{end}(:,1),split_arcs{end}(:,2),split_arcs{end}(:,3),'*');
        end
    end
    
    len_final_nodes = size(final_nodes,1);
    ones_nodes      = ones(len_final_nodes,4);
    for t=1:size(intersect_nodes,1)
        comp_nodes =  sum(ones_nodes.*intersect_nodes(t,:) - final_nodes,2);
        if min(abs(comp_nodes)) ~= 0
            final_nodes(end+1,:) = intersect_nodes(t,:);
            len_final_nodes = size(final_nodes,1);
            ones_nodes      = ones(len_final_nodes,4);
        end
    end
end


stacked_nodes = [start_nodes; end_nodes];
nodes = unique(stacked_nodes(:,1:3),'rows'); %Code only requires xyz coordinates

% Now double check
[who,where] = sort(nodes(:,3));
check_nodes = nodes(where,:);
diff_nodes = diff(check_nodes);
vessel_node_id = [];
figure; hold on;
for i=1:max(size(split_arcs))
    temp_start = split_arcs{i}(1,1:3);
    temp_end   = split_arcs{i}(end,1:3);
    where_start = find(start_nodes(:,1) == temp_start(1) & start_nodes(:,2) == temp_start(2) & start_nodes(:,3) == temp_start(3));
    where_end   = find(end_nodes(:,1) == temp_end(1) & end_nodes(:,2) == temp_end(2) & end_nodes(:,3) == temp_end(3));
    if isempty(where_start) || isempty(where_end)
        start_nodes(end+1,:) = split_arcs{i}(1,:);
        end_nodes(end+1,:)   = split_arcs{i}(end,:);
    end
    plot3(split_arcs{i}(:,1),split_arcs{i}(:,2),split_arcs{i}(:,3),'*')
end
hold off;

% NEED TO FIGURE OUT WHY THE END NODES DO NOT MATCH
if flag == 0
    nodes = [new_root(end,1:3); nodes];
    new_arcs = {};
    new_arcs{end+1} = new_root;
    for i=1:max(size(split_arcs))
        new_arcs{end+1} = split_arcs{i};
    end
    split_arcs = new_arcs;
end
end