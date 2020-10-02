% solves the model and plots solutions against data
% DRIVERBASIC.M
%
% This function takes in the data structure, and then runs the model until
% convergence has been met for all the model states.
function [rout,J,Lts,CV] = CVmodel(pars,data)

global ODE_TOL REL_TOL

pars = exp(pars);
Init = data.Init;
res_flag = data.res_flag;

% Resistances
Rs    = pars(1); % systemic periphery
Rp    = pars(2); % pulmonary periphery
Rava  = pars(3); % aortic valve
Rmva  = pars(4); % mitral valve
Rpva  = pars(5); % pulmonary valve
Rtva  = pars(6); % tricuspid valve
Rpv   = pars(7); % pulmonary veins
Rsv   = pars(8); % systemic veins

% Compliances
Csa = pars(9); % systemic arteries
Csv = pars(10); % systemic veins
Cpa = pars(11); % pulmonary arteries
Cpv = pars(12); % pulmonary veins

% Elatance Heart parameters
EMra   = pars(13); % max right atrial elastance
Emra   = pars(14); % min right atrial elastance
EMla   = pars(15); % max left ventricle elastance
Emla   = pars(16); % min right ventricle elastance
EMrv   = pars(17); % max left ventricle elatance
Emrv   = pars(18); % min left ventricle elatance
EMlv   = pars(19); % max left ventricle elatance
Emlv   = pars(20); % min left ventricle elatance

% Heart timing parameters
% Atrium
Trra = pars(21);      % end right atrium relaxation
tcra = Trra+pars(22); % right atrium begin contraction
Tcra = tcra+pars(23); % right atrium end contraction

Trla = Trra.*(1.01); % end left atrium relaxation MITCHEL SCALE!
tcla = tcra.*(1.05); % left atrium begin contraction
Tcla = Tcra; % left atrium end contraction

% Ventricle
Tcrv = pars(24);        % right ventricle contraction
Trrv = Tcrv + pars(25); % right ventricle relaxation

Tclv = Tcrv.*(0.95); % left ventricle contraction MITCHEL SCALE!
Trlv = Trrv; % left ventricle relaxation
    
% Initialize vectors for pressure solutions
CV.ppaS   = [];  
CV.pPAdS  = [];
ppaMS     = [];  
pPAMdS    = [];    
ppamS     = [];
pPAmdS    = [];
CV.praS   = [];
CV.pRAdS  = [];
praMS     = [];
pRAMdS    = [];
pramS     = [];
pRAmdS    = [];
CV.prvS   = [];
CV.pRVdS  = [];
prvMS     = [];
pRVMdS    = [];
prvmS     = [];
pRVmdS    = [];
CV.psaMS  = [];
CV.pSAMdS = [];
CV.psamS  = [];
CV.pSAmdS = [];
CV.ppvS   = [];
CV.pPWdS  = [];
CV.COsS   = [];
CV.COdS   = [];
CV.tdS    = [];

T  = data.T;
td = data.td;
dt = data.dt;
NC = data.NC; 

k1 = 1; % index of first time step in first period
k2 = round(T/dt)+k1; %index of last time step in first period
for i = 1:NC % go through loop NC times
        
    clear Vsa Vsv Vpa Vpv  Vra Vla Vrv Vlv
    clear psa psv ppa ppv  Era Era Ela pla Erv prv Elv plv
    clear qs qp  qmva qtva qsv qpv
        
    tdc = td(k1:k2);  %current time
        
    options=odeset('RelTol',REL_TOL, 'AbsTol',ODE_TOL);       %sets how accurate the ODE solver is
    sol = ode15s(@modelBasic,tdc,Init,options,pars,tdc(1),T); %solves the ODE uses modelBasic
    sols= deval(sol,tdc);                                     %interpolates solution at data times
        
    % extract solutions
    Vsa  = sols(1,:)';  % systemic arteries
    Vsv  = sols(2,:)';  % systemic veins
    Vpa  = sols(3,:)';  % pulmonary arteries
    Vpv  = sols(4,:)';  % pulmonary veins
    Vra  = sols(5,:)';  % right atrium
    Vla  = sols(6,:)';  % left atrium
    Vrv  = sols(7,:)';  % right ventricle
    Vlv  = sols(8,:)';  % left ventricle
        
    psa = Vsa/Csa; % systemic arteries
    psv = Vsv/Csv; % systemic veins
    ppa = Vpa/Cpa; % pulmonary arteries
    ppv = Vpv/Cpv; % pulmonary veins
        
    % Heart elastance
    for w = 1:length(tdc)
      Era(w)  = ElastanceAtrium(tdc(w)-tdc(1),EMra,Emra,Trra,tcra,Tcra,T); % Right atrium elastance
      Ela(w)  = ElastanceAtrium(tdc(w)-tdc(1),EMla,Emla,Trla,tcla,Tcla,T); % Left atrium elastance
            
      Erv(w)  = ElastanceVentricle(tdc(w)-tdc(1),EMrv,Emrv,Tcrv,Trrv);   % Right ventricle elastance
      Elv(w)  = ElastanceVentricle(tdc(w)-tdc(1),EMlv,Emlv,Tclv,Trlv);   % Right ventricle elastance
    end
    pra = Era'.*Vra; % Right atrium pressure
    pla = Ela'.*Vla; % Left atrium pressure
    prv = Erv'.*Vrv; % Right ventricle pressure
    plv = Elv'.*Vlv; % Left ventricle pressure
        
    % Flows defined by Ohm's Law
    qs    = (psa - psv)/Rs; % systemic periphery
    qp    = (ppa - ppv)/Rp; % pulmonary periiphery
        
    % Cardiac Output
    COs   = trapz(tdc, qs)/(tdc(end)-tdc(1));
    COp   = trapz(tdc, qp)/(tdc(end)-tdc(1));
        
    qmva = zeros(size(qs));
    if pla > plv % Mitral valve
      qmva = (pla-plv)/Rmva; % valve open
    end
        
    qtva = zeros(size(qs));
    if pra > prv  % Tricuspid valve
       qtva = (pra-prv)/Rtva; % valve open
    end
        
    % Venous flows
    qsv = zeros(size(qs));
    if psv > pra
       qsv = (psv-pra)/Rsv;
    end
        
    qpv = zeros(size(qs));
    if ppv > pla
       qpv = (ppv-pla)/Rpv;
    end
    
    % ppa
    CV.tdS    = [CV.tdS tdc(1:end-1)];
    CV.ppaS   = [CV.ppaS ppa(1:end-1)'];
    CV.pPAdS  = [CV.pPAdS data.pPA(1:end-1)'];
   
    ppaMS  = [ppaMS  max(ppa)];  
    pPAMdS = [pPAMdS max(data.pPA)];    
    ppamS  = [ppamS  min(ppa)];
    pPAmdS = [pPAmdS min(data.pPA)];
    
    % pra
    tcra_data = data.tcra+data.Trra;
    Tcra_data = tcra_data+data.Tcra;
    idmin = find( (tdc-tdc(1))>tcra_data,1)-7;
    idmax = find( (tdc-tdc(1))>Tcra_data,1)+7;
    
    CV.praS   = [CV.praS pra([1,idmin:idmax,end-1])'];
    CV.pRAdS  = [CV.pRAdS data.pRA([1,idmin:idmax,end-1])'];
    
    praMS  = [praMS  max(pra)];
    pRAMdS = [pRAMdS max(data.pRA)];
    pramS  = [pramS  min(pra)];
    pRAmdS = [pRAmdS min(data.pRA)];
    
    % prv
    CV.prvS   = [CV.prvS prv(1:end-1)'];
    CV.pRVdS  = [CV.pRVdS data.pRV(1:end-1)'];
    
    prvMS  = [prvMS  max(prv)];
    pRVMdS = [pRVMdS max(data.pRV)];
    prvmS  = [prvmS  min(prv)];
    pRVmdS = [pRVmdS min(data.pRV)];
    
    % psa
    CV.psaMS  = [CV.psaMS  max(psa)];
    CV.pSAMdS = [CV.pSAMdS data.pSAM];
    
    CV.psamS  = [CV.psamS  min(psa)];
    CV.pSAmdS = [CV.pSAmdS data.pSAm];
    
    % ppv
    CV.ppvS   = [CV.ppvS   mean(ppv)];
    CV.pPWdS  = [CV.pPWdS  min(data.pPW)];
    
    % CO
    CV.COsS   = [CV.COsS COs];
    CV.COdS   = [CV.COdS data.CO]; 
    
    Init = [Vsa(end) Vsv(end) Vpa(end) Vpv(end) Vra(end) Vla(end) Vrv(end) Vlv(end)];
    if i < NC
        k1 = k2; %sets last index of this loop as first index for next loop
        k2 = round(T/dt)+k1;
    end
end
% save solution
CV.tdS    = [CV.tdS tdc(end)]; 
CV.ppaS   = [CV.ppaS ppa(end)];
CV.pPAdS  = [CV.pPAdS data.pPA(end)];
CV.praS   = [CV.praS pra(end)];
CV.pRAdS  = [CV.pRAdS data.pRA(end)];
CV.prvS   = [CV.prvS prv(end)];
CV.pRVdS  = [CV.pRVdS data.pRV(end)];

M2   = round(NC/2);
N2   = round(length(CV.ppaS)/2);  % number of points per cycle x number of cycles / 2
NpRA = round(length(CV.praS)/2);

routS = [(ppaMS(M2:end)-pPAMdS(M2:end))./pPAMdS(M2:end)/sqrt(M2) ...
         (ppamS(M2:end)-pPAmdS(M2:end))./pPAmdS(M2:end)/sqrt(M2) ...
         (praMS(M2:end)-pRAMdS(M2:end))./pRAMdS(M2:end)/sqrt(M2) ...
         (pramS(M2:end)-pRAmdS(M2:end))./pRAmdS(M2:end)/sqrt(M2) ...
	     (prvMS(M2:end)-pRVMdS(M2:end))./pRVMdS(M2:end)/sqrt(M2) ...
         (prvmS(M2:end)-pRVmdS(M2:end))./pRVmdS(M2:end)/sqrt(M2) ...
         (CV.psaMS(M2:end)-CV.pSAMdS(M2:end))./CV.pSAMdS(M2:end)/sqrt(M2) ... 
         (CV.psamS(M2:end)-CV.pSAmdS(M2:end))./CV.pSAmdS(M2:end)/sqrt(M2) ...
         (CV.ppvS(M2:end) -CV.pPWdS(M2:end)) ./CV.pPWdS(M2:end)/sqrt(M2) ...
         (CV.COsS(M2:end) -CV.COdS(M2:end))  ./CV.COdS(M2:end)/sqrt(M2)];
LS   = 10;

if res_flag == 1
% Simulation 1 (Static values)   
     routT = [];
     LT    = 0;
elseif res_flag == 2
% Simulation 2 (Ventricular Pressure ONLY)
     routT = [(CV.prvS(N2:end)-CV.pRVdS(N2:end))./CV.pRVdS(N2:end)/sqrt(N2)];
     LT    = 1;
elseif res_flag == 3
% Simulation 3 (MPA Pressure ONLY)    
     routT = [(CV.ppaS(N2:end)-CV.pPAdS(N2:end))./CV.pPAdS(N2:end)/sqrt(N2)];
     LT    = 1;
elseif res_flag == 4
% Simulation 4 (RA Pressure ONLY) 
     routT = [(CV.praS(NpRA:end)-CV.pRAdS(NpRA:end))./CV.pRAdS(NpRA:end)/sqrt(NpRA)];
     LT    = 1;
elseif res_flag == 5
% Simulation 5 RV & MPA Time Varying      
     routT = [(CV.prvS(N2:end)-CV.pRVdS(N2:end))./CV.pRVdS(N2:end)/sqrt(N2) ...
              (CV.ppaS(N2:end)-CV.pPAdS(N2:end))./CV.pPAdS(N2:end)/sqrt(N2)];
     LT    = 2;
 elseif res_flag == 6    
% Simulation 6 (RA and RV Time-Varying)  
     routT   = [(CV.praS(NpRA:end)-CV.pRAdS(NpRA:end))./CV.pRAdS(NpRA:end)/sqrt(NpRA) ...
             (CV.prvS(N2:end)-CV.pRVdS(N2:end))./CV.pRVdS(N2:end)/sqrt(N2)];
     LT     = 2;
elseif res_flag == 7    
%Simulation 7 (RA and PA Time-Varying)  
    routT  = [(CV.praS(NpRA:end)-CV.pRAdS(NpRA:end))./CV.pRAdS(NpRA:end)/sqrt(NpRA) ...
              (CV.ppaS(N2:end)-CV.pPAdS(N2:end))./CV.pPAdS(N2:end)/sqrt(N2)];
    LT     = 2;
elseif res_flag == 8
% Simulation 8 ALL Time Varying        
    routT  = [(CV.prvS(N2:end)-CV.pRVdS(N2:end))./CV.pRVdS(N2:end)/sqrt(N2) ...
              (CV.ppaS(N2:end)-CV.pPAdS(N2:end))./CV.pPAdS(N2:end)/sqrt(N2) ...
              (CV.praS(NpRA:end)-CV.pRAdS(NpRA:end))./CV.pRAdS(NpRA:end)/sqrt(NpRA)];
    LT     = 3;
end
Lts  = [LT LS N2 M2];
L    = LT + LS; % number of components (sources of data)
rout = [routS routT]'/sqrt(L); 
J    = rout'*rout;

% figure(1);hold on;
% plot(tdc,data.pRV,'k',tdc,prv,'r');
% 
% figure(2);clf;
% plot(tdc,data.pRV,tdc,prv);
% pause;
end
