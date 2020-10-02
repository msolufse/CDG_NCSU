%
% call_flag: Defines which subroutines to call. Options include "model" for
% DriverBasic, "sens" for DriverBasic_sens, and "Stats" for RunStats.

load('Patient4.mat');

% you need to enter prompt with quotations 'Stats'
prompt = 'Which call_flag? sens, Stats: ';
call_flag = input(prompt);


switch call_flag  

    case 'Stats' 
    Stats = RunStats(data);
    info = 'Statistics for optimized parameters. from PostP4.';
    save('StatisticsP4_WSubsets.mat','Stats','info')

    k = 1;
    p = 1;
    for i = 1:1
        for j = 1:8
            R2_ppa(k)   = Stats{i,j}.R22(1);
            R2_prv(k)   = Stats{i,j}.R22(2);
            R2_pra(k)   = Stats{i,j}.R22(3);        
            J(k)        = Stats{i,j}.J;
            if Stats{i,j}.AIC == Inf
                continue
            else             
                AIC(k)  = Stats{i,j}.AIC;
                AICc(k) = Stats{i,j}.AICc;
                BIC(k)  = Stats{i,j}.BIC;
            end
            k = k+1;
        end 
    end 


    % Don't Run 'sens' part
    case 'sens'    
    % Sensitivity analysis
    Res = 8; 
    [pars, sens, Rsens, Isens]  = DriverBasic_sens(data,Res);   
    info = 'Sens results for Res 8';
    S = strcat('Sens',num2str(Res),'.mat');
    save (S,'sens','Rsens','Isens','pars','info');
 
        
end 
