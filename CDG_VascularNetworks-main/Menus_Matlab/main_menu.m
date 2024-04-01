clear all
close all

prompt         = {'Fontsize (PC: 12, Mac: 16)','Markersize (PC: 10, Mac: 14)','Linewidth (PC: 2, Mac: 3)'};
dlgtitle       = 'Select Figure Dimensions'; 
definput       = {'14','12','3'};
dims           = [1 65];
answer         = inputdlg(prompt,dlgtitle,dims,definput);

fontsize       = str2num(answer{1});
markersize     = str2num(answer{2});
linwidth       = str2num(answer{3});
cd 
if isempty(answer) == 0
    plotmarkers.fs = fontsize;
    plotmarkers.ms = markersize;
    plotmarkers.lwt= linwidth;
else
    return
end 

% Select an excel file (output from R)

prompt   = {'Select Data Set'};

dlgtitle = 'Dataset List';

d       = dir('ExcelData/*.xlsx');
liststr = {d.name};

prompt  = {'Select a file.'};
       
[indx,~] = listdlg('PromptString', prompt, 'SelectionMode','single', ...
    'ListString',liststr,'ListSize',[200,300],'Name','Patient Data Selection');

if isempty(indx) == 1
    return
end

patient_data = liststr(indx);


% Select either the Aorta or the Pulmonary Arteries to analyze
prompt_op = {'Select the system of interest.'};
call_flag_list = {'1. Aorta', ...
    '2. Pulmonary Arteries'};

[indx_op,tf] = listdlg('PromptString',prompt_op,...
    'SelectionMode','single','ListString',call_flag_list,...
    'ListSize',[200,300],'Name','Operation Selection',...
    'SelectionMode','single');


call_flag = call_flag_list{indx_op};
operations(call_flag, patient_data);

close all