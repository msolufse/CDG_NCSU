%% JUNCTION CORRECTION FOR BEFORE RADIUS DETECTION
% Last updated July 28, 2023 by Emma Slack

%% ABOUT JUNCTION CORRECTION CODE
% This program shifts the junction points in vessels to a more appropriate 
% location. We find the normalized distances between the first two
% daughters of each parent vessel and isolate the portion which is
% within 5% of each other (given that 100% is the maximum distance between the
% daughters). We call this 5% as the percent of togetherness. We then 
% average these sections in the daughters and append it to the parent. 
% We "smooth" the transitions from this new junction point to the daughters.
% Also, we account for trifurcations and vessels that are too short.

%% INPUT VALUE LOCATIONS - YOU SHOULD CHANGE THESE
% Percentages of togetherness and for the weighted average across the 5
% different lengths of vessels (Lines 195-208). You may also want to change
% the threshold value (Line 34).

function correct_junctions(patientData, connectivityData)

load(patientData);
load(connectivityData);

% Make a copy of vessel_details and the connectivity matrix so that you 
% still have the original data
connectivity2 = connectivity;
new_vessel_details = vessel_details;
new_conn_mat_num = conn;

%% CHANGE THRESHOLD VALUE
% The threshold variable is one less than the minimum number of data points
% we need in a vessel. If a vessel does not meet this minimum length, then
% we add it to its parent and reassign its daughters.
threshold = 10;

%% Moving Junctions
% We iterate backwards through our vessels as to elongate the vessels as
% much as possible.
for i = height(new_conn_mat_num):-1:1
    parent = new_conn_mat_num(i,1);
    daughter1 = new_conn_mat_num(i, 2);
    daughter2 = new_conn_mat_num(i, 3);
    daughter3 = new_conn_mat_num(i, 4);

    % Ignore junction correction if the vessel is terminal
    if daughter1 == 0
        continue;
    end

    % The xyz-coordinates of all three daughters
    daughter1_vals = flip(new_vessel_details{daughter1+2,2}(:,1:3));
    daughter2_vals = flip(new_vessel_details{daughter2+2,2}(:,1:3));
    daughter3_vals = flip(new_vessel_details{daughter3+2,2}(:,1:3));

    daughter1_radii = flip(new_vessel_details{daughter1+2,2}(:,4));
    daughter2_radii = flip(new_vessel_details{daughter2+2,2}(:,4));
    daughter3_radii = flip(new_vessel_details{daughter3+2,2}(:,4));

    % Given the length of the parent vessel, we will change the percent of
    % togetherness to account for this extra length. We also isolate the
    % case of the root vessel, as this vessel behaves different from the
    % rest.
    % if new_vessel_details{parent+2,3} <  0.25 * max_length
    %     percent_together = 0.10;
    %     percent_avgd = 0.20;
    % elseif new_vessel_details{parent+2,3} <  0.45 * max_length
    %     percent_together = 0.35;
    %     percent_avgd = 0.45;
    % elseif new_vessel_details{parent+2,3} <  0.65 * max_length
    %     percent_together = 0.45;
    %     percent_avgd = 0.65;
    % elseif parent == 0
    %     percent_together = 0.15;
    %     percent_avgd = 0.45;
    % else
    %     percent_together = 0.65;
    %     percent_avgd = 0.95;
    % end

    percent_together = 0.10;
    percent_avgd = 0.20;

    % If there exists a third daughter, then we will calculate the new
    % parent differently
    if daughter3 ~= 0
        % We override what the percent together should be to 20%
        percent_together = 0.2;

        % The differences cannot be larger than the minimum height of the
        % three daughters.
        num_diff = min([height(daughter1_vals), height(daughter2_vals),  height(daughter3_vals)]);

        % We find the differences between each pair of daughters
        diffs12 = daughter1_vals(1:num_diff, :) - daughter2_vals(1:num_diff, :);
        diffs13 = daughter1_vals(1:num_diff, :) - daughter3_vals(1:num_diff, :);
        diffs23 = daughter2_vals(1:num_diff, :) - daughter3_vals(1:num_diff, :);

        % For whichever difference is the least, we consider the average of
        % these two daughters and the remaining daughter for our
        % differences calculation.
        if diffs12(1,:) == min([diffs12(1,:), diffs13(1,:), diffs23(1,:)])
            conjoined_sib = (daughter1_vals(1:num_diff, :) + daughter2_vals(1:num_diff, :))/2;
            odd_sib = daughter3_vals(1:num_diff, :);
        elseif diffs13(1,:) == min([diffs12(1,:), diffs13(1,:), diffs23(1,:)])
            conjoined_sib = (daughter1_vals(1:num_diff, :) + daughter3_vals(1:num_diff, :))/2;
            odd_sib = daughter2_vals(1:num_diff, :);
        else
            conjoined_sib = (daughter2_vals(1:num_diff, :) + daughter3_vals(1:num_diff, :))/2;
            odd_sib = daughter1_vals(1:num_diff, :);
        end

        % Difference between the average of the two closest daughters and
        % the remaining daughter.
        diffs = conjoined_sib - odd_sib;

    else
        % We can only compare daughters of the same length, so we only look at
        % the daughters up to the length of the shorter daughter.
        num_diff = min(height(daughter1_vals), height(daughter2_vals));

        % Differrence calculation
        diffs = daughter1_vals(1:num_diff, :) - daughter2_vals(1:num_diff, :);
    end

    % Distance calculation
    dists = sqrt(sum(diffs.^2, 2));

    % If the daughters lie on top of each other, then we want to add all of 
    % this to the parent vessel.
    if isempty(find(dists, 1))
        new_parent = [flip(vessel_details{parent+2,2}(:,:)); flip(vessel_details{daughter1+2,2}(:,:))];
        new_vessel_details{parent+2,2} = flip(new_parent);
        continue;
    end

    % Vector of normalized distances between daughters
    norm_dists = dists/max(dists);

    % Vector of indices in norm_dists which have values less than 0.05
    % (i.e. which points are less than 5% away from their sibling). This
    % percentage can change, and should be changed for different lengths of
    % vessels.
    indices = find(norm_dists <= percent_together);

    if isempty(indices)
        indices = [height(num_diff)];
    end

    % If there is a third daughter, then we want to append the 
    if daughter3 ~= 0
        % We take the midpoint between three points (the centroid of the
        % triangle that is created between the three daughters). This is
        % equivalent to finding the average between the three daughters.
        addition = (daughter1_vals(1:indices(end), :) + daughter2_vals(1:indices(end), :) + daughter3_vals(1:indices(end), :))/3;
        new_end = addition(end, :);

        % Height of remaining sections in each daughter
        rad_height1 = height(daughter1_vals(indices(end)+1:end, :));
        rad_height2 = height(daughter2_vals(indices(end)+1:end, :));
        rad_height3 = height(daughter3_vals(indices(end)+1:end, :));

        % Vector of computed radius for height of remaining sections
        side1 = daughter1_radii(indices(end)+1:end, :);
        side2 = daughter2_radii(indices(end)+1:end, :);
        side3 = daughter3_radii(indices(end)+1:end, :);

        % Vector of the new daughter informatin
        new_radii = flip(new_vessel_details{parent+2,2}(end, 4));
        new_daughter1 = [new_end, new_radii; daughter1_vals(indices(end)+1:end, :), side1];
        new_daughter2 = [new_end, new_radii; daughter2_vals(indices(end)+1:end, :), side2];
        new_daughter3 = [new_end, new_radii; daughter3_vals(indices(end)+1:end, :), side3];

        % Reassigning the new daughter info into our table
        new_vessel_details{daughter1+2,2} = flip(new_daughter1);
        new_vessel_details{daughter2+2,2} = flip(new_daughter2);
        new_vessel_details{daughter3+2,2} = flip(new_daughter3);

        % Updating the parent vessel and adding it to our table
        parent_side = (daughter1_radii(1:indices(end),:)+daughter2_radii(1:indices(end), :)+daughter3_radii(1:indices(end), :))/3;
        new_parent = [flip(vessel_details{parent+2,2}(:,:)); addition, parent_side];
        new_vessel_details{parent+2,2} = flip(new_parent);
        continue;
    end

    % Appending the average of this section to the parent
    new_end = (daughter1_vals(indices(end), :) + daughter2_vals(indices(end), :))/2;
    addition = (daughter1_vals(1:indices(end), :) + daughter2_vals(1:indices(end), :))/2;
    parent_side = (daughter1_radii(1:indices(end),:)+daughter2_radii(1:indices(end), :))./2;
    new_parent = [flip(new_vessel_details{parent+2,2}(:,:)); addition, parent_side];
    new_vessel_details{parent+2,2} = flip(new_parent);

    % The "pretend parent" is the average between the two daughters for the
    % rest of the shortest daughter's length (this would be the parent if
    % we averaged 100% of the daughters).
    pretend_parent = (daughter1_vals(indices(end):num_diff, :) + daughter2_vals(indices(end):num_diff, :))/2;

    % Vector of indices from 5% to 10% of the max distance between
    % daughters
    indices_5to10 = find(norm_dists(indices(end)) <= percent_avgd) + indices(end);

    % If there are no points between 5% and 10% (or whatever our percentages
    % are), then we will not smooth.
    if isempty(indices_5to10)
        % Height of remaining sections in both daughters
        rad_height1 = height(daughter1_vals(indices(end):end, :));
        rad_height2 = height(daughter2_vals(indices(end):end, :));

        % Vector of computed radius for height of remaining sections
        side1 = daughter1_radii(1:rad_height1, :);
        side2 = daughter2_radii(1:rad_height2, :);

        % Reassigning the new daughter info into our table
        new_radii = new_vessel_details{parent+2,2}(end, 4);
        new_daughter1 = [daughter1_vals(indices(end):end, :), side1; new_end, new_radii];
        new_daughter2 = [daughter2_vals(indices(end):end, :), side2; new_end, new_radii];

    % If there are points within 5% and 10%, then we compute a weighted
    % average and "smooth" the daughters.
    else
        % We only want points after the portion reassigned to the parent in
        % our daughter vessels, so we calculate the indices as such.
        pretend_parent_index = indices_5to10 - indices_5to10(1)+1;

        % Weighted average for 5%-10% of the daughters
        new_d1_avg = [(daughter1_vals(indices_5to10(1):indices_5to10(end), :)*2 + pretend_parent(pretend_parent_index(1):pretend_parent_index(end), :))/3];
        new_d2_avg = [(daughter2_vals(indices_5to10(1):indices_5to10(end), :).*2 + pretend_parent(pretend_parent_index(1):pretend_parent_index(end), :))./3];

        % New xyz-coordinates for the daughters
        new_d1_xyz = [new_end; new_d1_avg; daughter1_vals(indices_5to10(end)+1:end, :)];
        new_d2_xyz = [new_end; new_d2_avg; daughter2_vals(indices_5to10(end)+1:end, :)];

        % Radii for the daughters
        flipped_d1 = flip(new_vessel_details{daughter1+2,2}(:, 4));
        flipped_d2 = flip(new_vessel_details{daughter2+2,2}(:, 4));

        % New daughter matrices with xyz-coordinates and radii
        new_daughter1 = [new_d1_xyz, flipped_d1(indices_5to10(1)-1:end)];
        new_daughter2 = [new_d2_xyz, flipped_d2(indices_5to10(1)-1:end)];
    end
    % Reassigning the new daughter information to our table
    new_vessel_details{daughter1+2,2} = flip(new_daughter1);
    new_vessel_details{daughter2+2,2} = flip(new_daughter2);
end

%% Updating Vessel Lengths
% Iterating through each vessel in our network
for h = 2:height(new_vessel_details)
    new_length = 0;

    % Iterating through the observation points
    for k = 1:height(new_vessel_details{h,2})-1

        % Finding the distance between two adjacent observation points
        diffs = new_vessel_details{h,2}(k,1:3) - new_vessel_details{h,2}(k+1,1:3);
        new_length = new_length + sqrt(sum(diffs.^2, 2));
    end

    % Updating the length in the new connectivity matrix and new_vessel_details
    new_vessel_details{h,3} = new_length;
end
% 
% conn_mat_small = [];
% 
% for i = 1:height(conn)
%     if conn(i, 1) ~= 0 & conn(i,2) ~= 0
%         conn_mat_small = [conn_mat_small; conn(i, :)];
%     end
% end
% for i = 1:height(conn_mat_small)
%     for j = 2:4
%         old_vessel_index = conn_mat_small(i, j);
%         if old_vessel_index ~= 0
%             new_vessel_index = find(conn_mat_small(:, 1) == old_vessel_index);
%             conn_mat_small(i, j) = new_vessel_index-1;
%         end
%     end
%     conn_mat_small(i, 1) = i-1;
% end
% connectivity = conn_mat_small;
% new_new_vessel_details = new_vessel_details;

vessel_details = new_vessel_details;
connectivity = connectivity2;

% Up to this point, I have been saving everything, but new_vessel_details
% and the new_conn_mat_num (the new connectivity matrix) are the most 
% important things to save.
save('junction_correct_b4R.mat', 'vessel_details', 'connectivity')

end