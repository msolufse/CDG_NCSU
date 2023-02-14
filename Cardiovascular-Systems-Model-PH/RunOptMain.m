%%% For Cluster %%%    
%pool = gcp('nocreate'); 
%parpool([2, 4]);

R = 2; % 1 static only & 2 static & dynamic

ind = {};
for ptnb = 1:9
     %INDMAP = [1:25]; % Residual 2 ALL PARAMETERS
     %IND    = [1:20];
     
     %INDMAP = [1 2 6 8 9:12 14 16 18 20 21:25]; % Residual 2 ALL PARAMETERS (same as before)
     %IND    = [1 2 6 8 9:12 14 16 18 20];
     %INDMAP = [1 2 6 8 9:12 14 16 18 20 22:25]; % no Trra (21)
     %IND    = [1 2 6 8 9:12 14 16 18 20]; 
     %INDMAP = [1 2 6   9:12 14 16 18 20 22:25]; % no Rsv (8) and Trra
     %IND    = [1 2 6   9:12 14 16 18 20]; 
     %INDMAP = [1 2 6   9:12 14    18 20 22:25]; % no Trra (21) no Rsv (8) no Emla (16)
     %IND    = [1 2 6   9:12 14    18 20]; 
     %INDMAP = [1 2 6   9:11 14    18 20 22:25]; % no Trra (21) no Rsv (8) no Emla (16), no Cpv (12)
     %IND    = [1 2 6   9:11 14    18 20];
     
     %INDMAP = [1 2 6 8 9:12 14    18 20 22:25]; % no Trra (21) no Emla (16)
     %IND    = [1 2 6 8 9:12 14    18 20]; 
     %INDMAP = [1 2 6 8 9:11 14    18 20 22:25]; % no Trra (21) no Emla (16) no Cpv (12)
     %IND    = [1 2 6 8 9:11 14    18 20]; 
     
     %INDMAP = [1 2 6 9:12 14 16 18 20 23:25]; % Residual 1 ALL PARAMETERS 
     %IND    = [1 2 6 9:12 14 16 18 20];
     
     %INDMAP = [1 2 6 9:12 14 16 18 20 23 25]; % Residual 1 no Tcrv
     %IND    = [1 2 6 9:12 14 16 18 20];
    
     %INDMAP = [1 2 6 9:12 14  18 20 23 25]; % Residual 1 no Tcrv (24) no Emla (16)
     %IND    = [1 2 6 9:12 14  18 20];
     
     INDMAP = [1 2 6 9:11 14  18 20 23 25]; % Residual 1 no Tcrv (24) no Emla (16) no Cpv (12)
     IND    = [1 2 6 9:11 14  18 20];
     
     % x0 = [Rs Rp Rava Rmva Rpva Rtva Rpv Rsv ...  % 1-8
     %      Csa Csv Cpa Cpv ...                     % 9-12
     %      EMra Emra EMla Emla ...                 % 13-16
     %      EMrv Emrv EMlv Emlv ...                 % 17-20
     %      Trra tcra Tcra ...                      % 21-23
     %      Tcrv Trrv]';                            % 24-25
          
     %task = 1;     % Run Sens
     %task = 2;     % Run Opt
     task  = 3;     % Plot Results
     %task = 4;     % Plot Sens
        
     par = 2; % Nom pars org (0), Opt pars (2)
     display(strcat('patient ',num2str(ptnb)));
     rng('Shuffle');
     rnum  = zeros(25,9);
     rt = -1+2*rand(length(IND),8);
%      rnum(IND,1:8) = rt;
     iternum = [1 2 3 4 5 6 7 8 9];
     for i = 1%1:9
         i
         iter = iternum(i);
         RunOpt_Clust(R,INDMAP,task,par,ptnb,iter,rnum(:,iter));  
     end;
    disp('end patient');
    
    clear IND INDMAP i iter iternum par rnum rt task
%     close all
end;

%delete(gcp('nocreate'));
