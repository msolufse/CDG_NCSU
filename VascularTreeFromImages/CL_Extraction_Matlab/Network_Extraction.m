%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Network_Extraction.m File 
% Author: M. J. Colebank
% Last update: 4/20/19
%
% Description: This code is provided the centerline file (CLFILE) obtained
% after using VMTK to obtain centerlines for a 3D geometry. The output of
% the file is as follows:
%
% alpha_final: The final scaling factor that is used to convert the voxel
% dimensions to mm (needed for mice Micro-CT images, but not typically
% needed for standard CT/MRI files).
%
% linfit, lin_r0: a linear model for the tapering of the arteries along the
% principal pathway: f(x) = lin_r0 - linfit*x
%
% expofit, expo_r0: an exponential model for the tapering of the arteries 
% along the principal pathway: f(x) = expo_r0*exp(-expofit*x)
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD [alpha_final,linfit,expofit,lin_r0,expo_r0] = Network_Extraction(CLFILE)
function [alpha_final] = Network_Extraction(CLFILE)
nodesname = strcat('CL','Nodes.mat');
dataname = strcat('CL','Data.mat');
%% New input to make sure we have the right information via .csv file
%Sometimes it has quotes, other times it doesn't....
radnamequotes = '"MaximumInscribedSphereRadius"';
%radnamequotes = 'Radius';
xnamequotes = '"Points:0"';
ynamequotes = '"Points:1"';
znamequotes = '"Points:2"';
radname = 'MaximumInscribedSphereRadius';
%radname = 'Radius';
xname = 'Points:0';
yname = 'Points:1';
zname = 'Points:2';
fid = fopen(CLFILE);
C = textscan(fid, '%s');
fclose(fid);
header = char(C{1}(1));
headernames = {};
headernameID = 1;
startID = 1;
for i=1:length(header)
    if strcmp(header(i),',')
        headernames{headernameID} = header(startID:i-1);
        headernameID = headernameID + 1;
        startID = i+1;
    end
end
headernames{headernameID} = header(startID:i);
%Now find the stuff you want
cl = csvread(CLFILE,1);
for j=1:headernameID
    if strcmp(headernames{j},radnamequotes) || strcmp(headernames{j},radname)
        r = cl(:,j);
    elseif strcmp(headernames{j},xnamequotes) || strcmp(headernames{j},xname)
        x = cl(:,j);
    elseif strcmp(headernames{j},ynamequotes) || strcmp(headernames{j},yname)
        y = cl(:,j);
    elseif strcmp(headernames{j},znamequotes) || strcmp(headernames{j},zname)
        z = cl(:,j);
    else
    end
end
cl = [r x y z];
%% User input to define if there needs to be scaling
% For this project: we want the scaling factor to be 0.86 (hardcoded into
% this version).
flagnodes = 0;
scalingfactor = 0;
prompt = 'Do you have the nodes for this Network? 1: Yes, 0: No\n';
flagnodes = input(prompt);
prompt = 'Does this network have a clamp? If so, input scaling below or put 0 for no scaling:\n';
scalingfactor = input(prompt);%0.86;
% if flagnodes > 0
%     scalingfactor = flagscale;
% else
%     scalingfactor = 0;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm 1: New function for finding unique points
[nodes,arcs,adjust] = find_nodes(cl); %Adjust can be used to set the root vessel to have z=0 at the inlet
vessel_keep = arcs;
save(nodesname,'nodes','arcs','adjust');%,'adjust','node_bif');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_vessels = length(vessel_keep);
% Plotting routine to check segments
figure(203)
hold on;
grid on;
names = [];
vessel_starts = zeros(length(vessel_keep),3);
vessel_ends = zeros(length(vessel_keep),3);
for i=1:length(nodes)-1
    if i<=length(vessel_keep)
        vessel_starts(i,:) = [vessel_keep{i}(1,1) vessel_keep{i}(1,2) vessel_keep{i}(1,3)];
        vessel_ends(i,:)   = [vessel_keep{i}(end,1) vessel_keep{i}(end,2) vessel_keep{i}(end,3)];
        plot3(vessel_keep{i}(:,1),vessel_keep{i}(:,2),vessel_keep{i}(:,3),'k');
    end
    plot3(nodes(i+1,1),nodes(i+1,2),nodes(i+1,3),'ro','MarkerSize',10,'MarkerFace','r');
end
plot3(nodes(1,1),nodes(1,2),nodes(1,3),'ro','MarkerSize',10,'MarkerFace','r');

%% Now match node points to their respective vessels
connectivity = cell(num_vessels,8); % This should probably just be dependent on the side
vessel_organized = cell(num_vessels,2); % of any n-furcations we have.. Will change in the future.
connectivity{1,1} = 'Parent Vessel';
connectivity{1,2} = 'Daughter 1';
connectivity{1,3} = 'Daughter 2';
connectivity{1,4} = 'Daughter 3';
connectivity{1,5} = 'Daughter 4';
connectivity{1,6} = 'Daughter 5';
connectivity{1,7} = 'Daughter 6';
connectivity{1,8} = 'Bifurcation IDs';
connectivity{1,9} = 'Location in Vessel Vector';
%% Algorithm 2: Define the main branch (i.e. main pulmonary artery)
% Just take the vessel that has the lowest point
for i=1:num_vessels
    if sum(vessel_keep{i}(1,:) == adjust)==4 || sum(vessel_keep{i}(end,:) == adjust)==4
       saveID = i;
       break;
    end
end

%% Start with a single Main pathway;
% Note: Letter indicates generation, number indicates branch
charID = 65;
VESSEL_NAME  = 0; % Will be used to name vessels START AT CURRNET SPOT
VESSEL_INDEX = 1;  % To increment the cell vessel_organized START AT CURRNET SPOT
CONN_INDEX   = 2;  % To define where we are in connectivity matrix
vessel_organized{1,1} = strcat(char(charID),num2str(VESSEL_NAME));
vessel_organized{1,2} = vessel_keep{saveID};
connectivity{CONN_INDEX,1} = vessel_organized{VESSEL_INDEX,1};
vessel_keep{saveID}   = [];
%% This is used to determine if your values are in pixels or in metric units
if scalingfactor == 0
    alpha_final=0;
elseif scalingfactor ~= 0
    alpha = find_scale(scalingfactor, vessel_organized{1,2});
    alpha_final = alpha;
    nodes = alpha.*nodes;
    vessel_starts = vessel_starts.*alpha;
    for i=1:length(vessel_keep)
        vessel_keep{i} = alpha.*vessel_keep{i};
    end
    vessel_organized{1,2} = alpha.*vessel_organized{1,2};
end

%% Algorithm 3: We will now find the the daughters
% start_node = nodes(1,:);
start_vessel = vessel_organized{1,2};
charID = charID+1;
newlengthflag = 1;
[vessel_organized, ~, connectivity, ~, ~, ~, ~] =  ...
    find_network_cl(start_vessel,vessel_keep, vessel_organized, connectivity, VESSEL_NAME, VESSEL_INDEX, ...
    CONN_INDEX, charID, vessel_starts);

%%

%% Get rid of empty components
N = length(vessel_organized);
nonempty_vessel_organized = {}; nonempty_connectivity = {};
nonempty_connectivity{1,1} = 'Parent Vessel';
nonempty_connectivity{1,2} = 'Daughter 1';
nonempty_connectivity{1,3} = 'Daughter 2';
nonempty_connectivity{1,4} = 'Daughter 3';
nonempty_connectivity{1,5} = 'Daughter 4';
nonempty_connectivity{1,6} = 'Daughter 5';
nonempty_connectivity{1,7} = 'Daughter 6';
nonempty_connectivity{1,8} = 'Bifurcation IDs';
nonempty_connectivity{1,9} = 'Location in Vessel Vector';
for i=1:N %a check to ensure there are no empty entries
   if ~isempty(vessel_organized{i})
      nonempty_vessel_organized{end+1,1} = vessel_organized{i,1};
      nonempty_vessel_organized{end,2}   = vessel_organized{i,2};
   end
   if i < N
       if ~isempty(connectivity{i+1,1})
          nonempty_connectivity{end+1,1}     = connectivity{i+1,1};
          nonempty_connectivity{end,2}       = connectivity{i+1,2};
          nonempty_connectivity{end,3}       = connectivity{i+1,3};
          nonempty_connectivity{end,4}       = connectivity{i+1,4};
          nonempty_connectivity{end,5}       = connectivity{i+1,5};
          nonempty_connectivity{end,6}       = connectivity{i+1,6};
          nonempty_connectivity{end,7}       = connectivity{i+1,7};
          nonempty_connectivity{end,8}       = connectivity{i+1,8};
          nonempty_connectivity{end,9}       = connectivity{i+1,9};
       end
   end
end
% Add one more loop because size of connectivity is greater than that of
% vessel_organized
if ~isempty(connectivity{end,1}) && ...
        ~strcmp(nonempty_connectivity{end,1},connectivity{end,1})
          nonempty_connectivity{end+1,1}     = connectivity{end,1};
          nonempty_connectivity{end,2}       = connectivity{end,2};
          nonempty_connectivity{end,3}       = connectivity{end,3};
          nonempty_connectivity{end,4}       = connectivity{end,4};
          nonempty_connectivity{end,5}       = connectivity{end,5};
          nonempty_connectivity{end,6}       = connectivity{end,6};
          nonempty_connectivity{end,7}       = connectivity{end,7};
          nonempty_connectivity{end,8}       = connectivity{end,8};
          nonempty_connectivity{end,9}       = connectivity{end,9};
end

%% ADD A LOOP HERE TO MAKE TRIFURCATIONS INTO BIFURCATIONS
% WONT NEED IN FUTURE VERSIONS


%% Now reassign the two entities
vessel_organized = nonempty_vessel_organized;
connectivity = nonempty_connectivity;
N = length(vessel_organized);
figure(2);clf;hold on;grid on;
for i=1:size(vessel_organized,1)
     plot3(vessel_organized{i,2}(:,1),vessel_organized{i,2}(:,2),vessel_organized{i,2}(:,3))
    text(vessel_organized{i,2}(end,1),vessel_organized{i,2}(end,2),vessel_organized{i,2}(end,3),vessel_organized{i,1},'FontSize',20)
end
hold off;
%% Append length, radii, ect. to a struct.
vessel_details = cell(length(vessel_organized),10);
vessel_details{1,1} = 'NAME';
vessel_details{1,2} = 'CL FILE';
vessel_details{1,3} = 'LENGTH [mm]';
vessel_details{1,4} = 'MEAN RADII [mm]';
vessel_details{1,5} = 'MEDIAN RADII';
vessel_details{1,6} = 'TAPER';
vessel_details{1,7} = 'INLET [mm]';
vessel_details{1,8} = 'OUTLET [mm]';
vessel_details{1,9} = 'R_STD';
vessel_details{1,10} = 'PARENT';
vessel_details{1,11} = 'DAUGHTERS';
vessel_details{2,1} = vessel_organized{1,1};
vessel_details{2,2} = vessel_organized{1,2};
%% Algorithm 4: Update values of interest (except length) AND RADII SHOULD NOT GO HERE
% [vessel_details,alpha] = update_details(vessel_details, 2, N, connectivity, 0);
[vessel_details] = update_details_cl(vessel_details, 2, N, connectivity);
for i=3:length(vessel_organized)+1
    vessel_details{i,1} = vessel_organized{i-1,1};
    vessel_details{i,2} = vessel_organized{i-1,2};
%     [vessel_details,alpha] = update_details_cl(vessel_details,i, N, connectivity, alpha);
    [vessel_details] = update_details_cl(vessel_details,i, N, connectivity);
end


% %% Calculate centroids and use them to extend length measurements
% jumps_morethan2 = diff(trif);
% where_morethan2 = find(jumps_morethan2 > 2);%figure out how many trifurcations we have
% if isempty(where_morethan2) && ~isempty(trif)
%     where_morethan2 = 1;
% end
% jumps_only2 = diff(false_bif);
% where_only2 = find(jumps_only2 > 2);%figure out how many trifurcations we have
% if isempty(where_only2) && ~isempty(false_bif)
%     where_only2 = 1;
% end

%% A check to see if anything is repeated:


figure(202)
set(gca,'FontSize', 32)
display(connectivity)
display(vessel_details)
save(dataname,'connectivity','vessel_details');

figure(399);clf;
hold on;
set(gca,'FontSize',32);
for i=2:size(vessel_details,1)
    if strcmp(vessel_details{i,1},'EMPTY VESSEL')
    else
        m = size(vessel_details{i,2},1);
        half = round(m/2)+1;
        plot3(vessel_details{i,2}(:,1),vessel_details{i,2}(:,2),(-1).*vessel_details{i,2}(:,3),'LineWidth',2)
%         text(vessel_details{i,2}(half,1),vessel_details{i,2}(half,2),(-1).*vessel_details{i,2}(half,3), vessel_details{i,1},'FontSize',30)
    end
end

%% Algorithm 6: Find the principal Pathway and calculate tapering
% [principal_pathway1, principal_pathway2] = find_principal_pathway(vessel_details);
% principal_tot = [principal_pathway1 principal_pathway2];
% 
% % To find left or right, we compare the first component of the centerline
% % file. More negative -> right branch.
% if vessel_details{principal_pathway1(2),2}(end,1) > vessel_details{principal_pathway2(2),2}(end,1)
%     leftbranch = principal_pathway1;
%     rightbranch = principal_pathway2;
% else
%     leftbranch = principal_pathway2;
%     rightbranch = principal_pathway1;
% end
% left_cum_length = 0;
% left_radii = [];
% right_cum_length = 0;
% right_radii = [];
% left_x = 0;
% right_x = 0;
% %%
% for i=1:length(leftbranch)
%    figure(399);
%    hold on;
%    id = leftbranch(i);
%    clpath = vessel_details{id,2}(end:-1:1,:);
%    m = size(clpath,1);
%    half = round(m/2)+1;
%    plot3(clpath(:,1),clpath(:,2),(-1).*clpath(:,3),'ok','MarkerSize',5,'MarkerFace','k')
%    text(clpath(half,1),clpath(half,2),(-1).*clpath(half,3), vessel_details{id,1},'FontSize',30)
%    for k=1:size(clpath,1)-1
%       left_cum_length = left_cum_length + sqrt( (clpath(k,1) - clpath(k+1,1)).^2 ...
%          +(clpath(k,2) - clpath(k+1,2)).^2 + (clpath(k,3) - clpath(k+1,3)).^2);
%      left_radii(end+1) = clpath(k,4);
%    end
%    left_radii(end+1) = clpath(end,4);
% %    xspace = left_x:left_x+size(clpath,1)-1;
%    left_x = left_x+size(clpath,1);
% end
% M = length(left_radii)-1;
% xspace = 0:left_cum_length/M:left_cum_length;
% figure(10000);
% hold on;
% title('LEFT PATHWAY')
% plot(xspace,left_radii,'o');
% hold off;
% figure(10000);
% hold on;
% %Fit curves using the entire spatial sequence of radii
% [r0_left_lin,fitted_linear_left] = find_taper_linear(left_radii,left_cum_length);
% [r0_left_expo, fitted_expo_left] = find_taper_expo(left_radii,left_cum_length);
% left_linear_fit = @(x) r0_left_lin - fitted_linear_left.*x;
% left_expo_fit   = @(x) r0_left_expo.*exp(-fitted_expo_left.*x);
% h(1) = plot(xspace,left_linear_fit(xspace),'r','LineWidth',3);
% h(2) = plot(xspace,left_expo_fit(xspace),'--r','LineWidth',3);
% legend(h,{'Lin fitted','Expo fitted'})
% hold off;
% %%
% for i=1:length(rightbranch)
%    figure(399);
%    hold on;
%    id = rightbranch(i);
%    clpath = vessel_details{id,2}(end:-1:1,:);
%    m = size(clpath,1);
%    half = round(m/2)+1;
%    plot3(clpath(:,1),clpath(:,2),(-1).*clpath(:,3),'ok','MarkerSize',5,'MarkerFace','k')
%    text(clpath(half,1),clpath(half,2),(-1).*clpath(half,3), vessel_details{id,1},'FontSize',30)
%    for k=1:size(clpath,1)-1
%       right_cum_length = left_cum_length + sqrt( (clpath(k,1) - clpath(k+1,1)).^2 ...
%          +(clpath(k,2) - clpath(k+1,2)).^2 + (clpath(k,3) - clpath(k+1,3)).^2);
%      right_radii(end+1) = clpath(k,4);
%    end
%    right_radii(end+1) = clpath(end,4);
%    right_x = right_x+size(clpath,1);
% end
% M = length(right_radii)-1;
% xspace = 0:right_cum_length/M:right_cum_length;
% figure(10001);
% hold on;
% title('RIGHT PATHWAY')
% plot(xspace,right_radii,'o');
% hold off;
% right_start = mean(vessel_details{rightbranch(2),2}(:,4));
% right_end = mean(vessel_details{rightbranch(end),2}(:,4));
% figure(10001);
% hold on;
% %Fit curves using the entire spatial sequence of radii
% [r0_right_lin,fitted_linear_right] = find_taper_linear(right_radii,right_cum_length);
% [r0_right_expo, fitted_expo_right] = find_taper_expo(right_radii,right_cum_length);
% right_linear_fit = @(x) r0_right_lin - fitted_linear_right.*x;
% right_expo_fit   = @(x) r0_right_expo.*exp(-fitted_expo_right.*x);
% h(1) = plot(xspace,right_linear_fit(xspace),'r','LineWidth',3);
% h(2) = plot(xspace,right_expo_fit(xspace),'--r','LineWidth',3);
% legend(h,{'Lin fitted','Expo fitted'})
% hold off;
% linfit  = [fitted_linear_left; fitted_linear_right];
% expofit = [fitted_expo_left; fitted_expo_right];
% lin_r0 = [r0_left_lin; r0_right_lin];
% expo_r0 = [r0_left_expo; r0_right_expo];
% 


end

