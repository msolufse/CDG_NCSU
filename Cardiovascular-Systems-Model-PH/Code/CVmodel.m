% solves the model and plots solutions against data
% DRIVERBASIC.M
%
% This function takes in the data structure, and then runs the model until
% convergence has been met for all the model states.
function [out_P, out_V, out_q] = CVmodel(pars,data)

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


    qava = zeros(size(qs));
    ids = 1:length(qs);
    ids = ids(plv>psa); 
    for j=ids % flow through aortic valve
        qava(j) = (plv(j)-psa(j))/Rava;
    end
    
    qpva = zeros(size(qs));
    ids = 1:length(qs);
    ids = ids(prv>ppa);
    for j=ids % flow through pulmonary valve
        qpva(j) = (prv(j)-ppa(j))/Rpva;
    end
        
    qmva = zeros(size(qs));
    ids = 1:length(qs);
    ids = ids(pla>plv);
    for j=ids % Mitral valve
      qmva = (pla-plv)/Rmva; % valve open
    end
        
    qtva = zeros(size(qs));
    ids = 1:length(qs);
    ids = ids(pra>prv);
    for j=ids  % Tricuspid valve
       qtva = (pra-prv)/Rtva; % valve open
    end
        
    
    qsv = zeros(size(qs));
    ids = 1:length(qs);
    ids = ids(psv>pra);
    for j=ids % Systemic venous flow
       qsv = (psv-pra)/Rsv;
    end
        
    qpv = zeros(size(qs));
    ids = 1:length(qs);
    ids = ids(ppv>pla);
    for j=ids % Pulmonary venous flow
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


out_P = [pla plv psa psv pra prv ppa ppv];
out_V = [Vla Vlv Vsa Vsv Vra Vrv Vpa Vpv];
out_q = [qmva qava qs qsv qtva qpva qp qpv];
end
