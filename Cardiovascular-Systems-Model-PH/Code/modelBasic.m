function xdot = modelBasic(t,y,pars,tst,T)

Vsa    = y(1);  % systemic artery volume
Vsv    = y(2);  % systemic venous volume
Vpa    = y(3);  % pulmonary artery volume
Vpv    = y(4);  % pulmonary venous volume
Vra    = y(5);  % right atria volume
Vla    = y(6);  % left atria volume
Vrv    = y(7);  % right ventricular volume
Vlv    = y(8);  % left ventricular volume

%% PARAMETERS
% Resistances
Rs   = pars(1);  % systemic arteries
Rp   = pars(2);  % pulmonary arteries
Rava = pars(3);  % aortic valve =
Rmva = pars(4);  % mitral valve
Rpva = pars(5);  % pulmonary valve
Rtva = pars(6);  % tricuspid valve
Rpv  = pars(7);  % pulmonary veins
Rsv  = pars(8);  % systemic veins

% Compliances
Csa = pars(9); % systemic arteries
Csv	= pars(10); % systemic veins
Cpa = pars(11); % pulmonary arteries
Cpv = pars(12); % pulmonary veins

% Heart
% Elastance
EMra = pars(13); % right atrium max 
Emra = pars(14); % right atrium min 
EMla = pars(15); % left atrium max 
Emla = pars(16); % left atrium min 
EMrv = pars(17); % right ventricle max 
Emrv = pars(18); % right ventricle min 
EMlv = pars(19); % left ventricle max
Emlv = pars(20); % left ventricle min 

% Timing
Trra = pars(21); % end right atrium relaxation
tcra = Trra+pars(22); % right atrium begin contraction
Tcra = tcra+pars(23); % right atrium end contraction

Trla = Trra.*(1.01);  % end left atrium relaxation
tcla = tcra.*(1.05);   % left atrium begin contraction
Tcla = Tcra;  % left atrium end contraction

Tcrv = pars(24); % right ventricle contraction
Trrv = Tcrv+pars(25); % right ventricle relaxation
Tclv = Tcrv.*(0.95);     % left ventricle contraction
Trlv = Trrv;     % left ventricle relaxation

% Timevarying elastance
Era = ElastanceAtrium(t-tst,EMra,Emra,Trra,tcra,Tcra,T); % Elatance right atrium
Ela = ElastanceAtrium(t-tst,EMla,Emla,Trla,tcla,Tcla,T); % Elastance left atrium
Erv = ElastanceVentricle(t-tst,EMrv,Emrv,Tcrv,Trrv);     % Elastance right ventricle
Elv = ElastanceVentricle(t-tst,EMlv,Emlv,Tclv,Trlv);     % Elastance left ventricle

psa = Vsa/Csa; % systemic arteries 
psv = Vsv/Csv; % systemic veins 
ppa = Vpa/Cpa; % pressure pulmonary artery
ppv = Vpv/Cpv; % pressure pulmonary vein

% Pressure 
pra  = Era*Vra; % pressure right atria
prv  = Erv*Vrv; % pressure right ventricle
pla  = Ela*Vla; % pressure left atria
plv  = Elv*Vlv; % pressure left ventricle

if psv > pra
    qsv= (psv-pra)/Rsv;
else
  qsv = 0;
end;

% Flow through valves
if pra > prv  % flow through tricuspid valve
    qtva = (pra-prv)/Rtva; % valve open
else
    qtva = 0;             % valve closed
end

if pla > plv % flow through mitral valve
    qmva = (pla-plv)/Rmva; % valve open
else
    qmva = 0;             % valve closed
end

if plv > psa % flow through aortic valve 
    %dqava = (plv-psa - qava*Rava)/Lava; % valve open
    qava = (plv-psa)/Rava;
else
%   dqava = 0;                       % valve closed
    qava  = 0;
end

if prv > ppa % flow through pulmonary valve
%   dqpva = (prv-ppa-qpva*Rpva)/Lpva; % valve open
    qpva = (prv-ppa)/Rpva;
else
%   dqpva = 0;                     % valve closed
    qpva  = 0;
end

% Peripheral flows
qs   = (psa - psv) / Rs; % Systemic periphery
qp   = (ppa - ppv) / Rp; % Pulmonary periphery

% Venous flows (with venous valves - flow cannot go backwards into the veins)
if ppv > pla
  qpv = (ppv-pla)/Rpv;
else
  qpv = 0;
end;

if psv > pra
    qsv= (psv-pra)/Rsv;
else
  qsv = 0;
end;

% Differential Equations
dVsa = qava - qs;
dVsv = qs   - qsv;
dVpa = qpva - qp;
dVpv = qp   - qpv;
dVra = qsv  - qtva; 
dVla = qpv  - qmva;
dVrv = qtva - qpva;
dVlv = qmva - qava;

xdot = [dVsa;...
        dVsv;...
        dVpa;...
        dVpv;...
        dVra;...
        dVla;...
        dVrv;...
        dVlv;];
end