
function Call_Model(PATIENT,RESIDUAL) 
       
% If you had patient data, you would call "load_global_PH"
% [pars,Vunstr,Init,low,hi,data] = load_global_PH(data,ptnb);


% Here, we load the optimal parameters and patient data files from the
% manuscript.
% If you want to plot the NOMINAL predictions, load "pars" instead of
% "optpars"
if PATIENT<10
fname = strcat('OUT/P',num2str(PATIENT),'_R',num2str(RESIDUAL),'_Final.mat');
load(fname,'data','optpars');
elseif PATIENT==10
    [optpars,Vunstr,Init,data] = load_global_CTR;
    data.res_flag = RESIDUAL;
else
    error('Patient file does not exist.')
end


% Parameters are saved in log-scale, so we need to exponetiate

DriverBasic_plot(optpars,data,PATIENT);

end
