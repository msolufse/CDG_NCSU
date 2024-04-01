close all; clear all;
%% NOTES
% Things to change in this code:
% Data being loaded in (line 8)
% Name of the excel file (line 116 & 128)
% Name of the new data set being save post changepoint analysis (line 159)
%% load data
data     = load('CLData.mat');
vessel_details  = data.vessel_vessel_details;

%% Understand how data is extracted using the first vessel

% Radius   = flip(vessel_details{2,2}(:,4)./10);
% Length   = [linspace(0,vessel_details{2,3}./10,length(Radius))]';
% x        = flip(vessel_details{2,2}(:,1));
% y        = flip(vessel_details{2,2}(:,2));
% z        = flip(vessel_details{2,2}(:,3));
% Parent   = vessel_details{2,10};
% Daughter = vessel_details{2,11};
% Sheet    = vessel_details{2,1};
% % names    = ["Radius" "Length" "x" "y" "z" "Parent" "Daughter"];
% 
% table    = table(Radius, Length, x, y, z);
% % table.Parent = repmat({''},length(table.Radius),1);
% table(1,'Parent') = {Parent};
% % table.Daughter = repmat({''},length(table.Radius),1);
% table(1,'Daughter') = {Daughter};
% %% Write to excel sheet
% 
% writetable(table,'TestFile.xlsx','Sheet',3);

%% automate for entire dataset
% extract the correct information from each vessel (length, radius, xyz
% coordinates, parent & daughter(s)). Note that we have to flip the radius
% and xyz components because VMTK works backwards from the outlets to the
% inlet(s).

for i = 2:length(vessel_details(:,1))
    Radius       = flip(vessel_details{i,2}(:,4)./10); %converted from mm to cm
    Length       = [linspace(0,vessel_details{i,3}./10,length(Radius))]'; %converted from mm to cm
    x            = flip(vessel_details{i,2}(:,1));
    y            = flip(vessel_details{i,2}(:,2));
    z            = flip(vessel_details{i,2}(:,3));
    Parent       = vessel_details{i,10};
    Daughter     = vessel_details{i,11};
    Sheet_Name   = vessel_details{i,1}; %Takes the name of the vessel and makes it the corresponding sheet name in excel
    comma        = strfind(Daughter,','); %detects where commas are in the string defining the daughters
    terminal     = strfind(Daughter,'TERMINAL'); %detects which vessels have terminal under the daughter column 
    blank_dghtr  = 'NONE';
    T            = table(Radius, Length, x, y, z); %table with all numerical values

% Now make sure that each daughter gets put into their own column so that
% both daughters can be analyzed in the R changepoint code. This code also
% considers trifurcations. If a vessel is terminal it will have a column
% labeled "Daughter1" and the first cell will have 'TERMINAL' as the entry

    if comma(1)==3 && comma(2)==6 && comma(3)==7
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:5);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {blank_dghtr};
    elseif comma(1)==3 && comma(2)==8 && comma(3)==9
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:7);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {blank_dghtr};
    elseif comma(1)==3 && comma(2)==6 && comma(3)==9
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:5);
        Daughter3 = Daughter(7:8);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {Daughter3};
    elseif comma(1)==3 && comma(2)==6 && comma(3)==10
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:5);
        Daughter3 = Daughter(7:9);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {Daughter3};
    elseif comma(1)==3 && comma(2)==6 && comma(3)==11
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:5);
        Daughter3 = Daughter(7:9);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {Daughter3};
    elseif comma(1)==3 && comma(2)==6 && comma(3)==12
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:5);
        Daughter3 = Daughter(7:10);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {Daughter3};
    elseif comma(1)==3 && comma(2)==7 && comma(3)==8
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:6);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {blank_dghtr};
    elseif comma(1)==3 && comma(2)==7 && comma(3)==11
        Daughter1 = Daughter(1:2);
        Daughter2 = Daughter(4:6);
        Daughter3 = Daughter(8:10);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {Daughter3};
    elseif comma(1)==4 && comma(2)==8 && comma(3)==9
        Daughter1 = Daughter(1:3);
        Daughter2 = Daughter(5:7);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {blank_dghtr};
    elseif comma(1)==4 && comma(2)==9 && comma(3)==10
        Daughter1 = Daughter(1:3);
        Daughter2 = Daughter(5:8);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {blank_dghtr};
    elseif comma(1)==4 && comma(2)==8 && comma(3)==12
        Daughter1 = Daughter(1:3);
        Daughter2 = Daughter(5:7);
        Daughter3 = Daughter(9:11);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {Daughter3};
    elseif comma(1)==5 && comma(2)==10 && comma(3)==11
        Daughter1 = Daughter(1:4);
        Daughter2 = Daughter(6:9);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {blank_dghtr};
    elseif comma(1)==5 && comma(2)==10 && comma(3)==15
        Daughter1 = Daughter(1:4);
        Daughter2 = Daughter(6:9);
        Daughter3 = Daughter(11:14);
        T(1,'Daughter1') = {Daughter1};
        T(1,'Daughter2') = {Daughter2};
        T(1,'Daughter3') = {Daughter3};
    elseif terminal(1)==1
        T(1,'Daughter1') = {Daughter};
        T(1,'Daughter2') = {blank_dghtr};
        T(1,'Daughter3') = {blank_dghtr};
    end

% Create column with the parent vessel name
    T(1,'Parent') = {Parent};
% Export to excel where each vessel has it's own spreadsheet. The columns
% are in the following order: Radius, Length, xyz, Daughter 1, 2, etc.,
% Parent. In R two additional columns will be added that are the location
% of the changepoints and the slopes of the piecewise functions,
% respectively. The excel file is written to the current working directory
% unless otherwise specified.
    writetable(T,'Healthy_patient_data.xlsx', 'Sheet',Sheet_Name);
end 

%% Export back to matlab vessel_vessel_details structure
% write an automated loop that will overwrite the radius, length, and xyz
% coordinates that were redefined in R. We need to do this because
% this is what the create_fluids code will use to define the connectivity
% for large networks that cannot be defined by hand.


% Read in Excel file and make into cell
fn          = 'Pluto_postostium.xlsx';        
tBC         = [];
opt         = detectImportOptions(fn);
shts        = sheetnames(fn);
new_vessel_details = cell(length(shts),2);

for j = 1:length(shts)
    new_vessel_details{j,1} = [shts(j)];
    new_vessel_details{j,2} = [tBC;readtable(fn,opt,'Sheet',shts(j))];
end

vessel_details{1,12} = [];
vessel_details{1,12} = 'CP COORDINATES';
% vessel_details{2,12} = zeros(nnz(~isnan(new_vessel_details{s,2}.cpx)),3);


% Overwrite old vessel_vessel_details with new components
for s = 1:height(new_vessel_details)
    vessel_details{s+1,12}      = zeros(1,3);
    vessel_details{s+1,2}       = zeros(length(new_vessel_details{s,2}.Radius),4);
    vessel_details{s+1,2}(:,4)  = new_vessel_details{s,2}.Radius; %i+1 because vessel_details has headers and new_vessel_details does not
    vessel_details{s+1,2}(:,1)  = new_vessel_details{s,2}.x;
    vessel_details{s+1,2}(:,2)  = new_vessel_details{s,2}.y;
    vessel_details{s+1,2}(:,3)  = new_vessel_details{s,2}.z;
    vessel_details{s+1,3}       = new_vessel_details{s,2}.Length(end);
    vessel_details{s+1,12}(1,1) = str2double(new_vessel_details{s,2}.cpx(1));
    vessel_details{s+1,12}(2,1) = str2double(new_vessel_details{s,2}.cpx(2));
    vessel_details{s+1,12}(1,2) = str2double(new_vessel_details{s,2}.cpy(1));
    vessel_details{s+1,12}(2,2) = str2double(new_vessel_details{s,2}.cpy(2));
    vessel_details{s+1,12}(1,3) = str2double(new_vessel_details{s,2}.cpz(1));
    vessel_details{s+1,12}(2,3) = str2double(new_vessel_details{s,2}.cpz(2));
%     vessel_details{s+1,12}(:,1) =
%     new_vessel_details{s,2}.cpx(1:nnz(~isnan(new_vessel_details{s,2}.cpx)));
%     %presostium
%     vessel_details{s+1,12}(:,2) =
%     new_vessel_details{s,2}.cpy(1:nnz(~isnan(new_vessel_details{s,2}.cpy)));
%     %presostium
%     vessel_details{s+1,12}(:,3) =
%     new_vessel_details{s,2}.cpz(1:nnz(~isnan(new_vessel_details{s,2}.cpz))); %preostium
    rtop = str2double(new_vessel_details{s,2}(1,15));
    if isnan(rtop) == true
        vessel_details{s+1,7} = (mean(new_vessel_details{s,2}.Radius))*10;
        vessel_details{s+1,8} = (mean(new_vessel_details{s,2}.Radius))*10;
    else
        vessel_details{s+1,7}       = new_vessel_details{s,2}.rtop(1)*10; % comment out if pre ostium removal
        vessel_details{s+1,8}       = new_vessel_details{s,2}.rbot(1)*10; % comment out if pre ostium removal
    end 
end

vessel_vessel_details = vessel_details;
connectivity   = data.connectivity; %uncomment if you have loaded the
% data in the matlab file not imported it
save('Orion_postostium_delete.mat','connectivity','vessel_vessel_details');

% now you can use the new .mat file created in line 159 to run 
% create_data_fluids.m to create the connectivity and dimension 
% files for running the fluids code