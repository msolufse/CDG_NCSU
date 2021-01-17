 % solves the model and plots solutions against data
% DRIVERBASIC.M

function Driver_treatment(PATIENT,TREATMENT)%data,xopt,I,fig
global ODE_TOL REL_TOL

filename = strcat('Pat',num2str(PATIENT),'_res2_opt.mat');
%% Load in the data and parameters to be used in the model.
% load the optimized file
load(filename,'data','pars');
Init = data.Init;
pars = exp(pars); %parameters are the optimize values mix in with non optimize parameters

%% Run the specific treatment for each patients
INDMAP = [1 2 9 11]; %% Rs, Rp, Csa, and Cpa parameters
switch TREATMENT
    % PAH & CTEPH treatment - vasodilators
    case 1
        par_change = [0.7 0.6 1.3 1.4];
    case 2
        par_change = [0.8 0.6 1.2 1.4];
    case 3
        par_change = [0.9 0.6 1.1 1.4];
    case 4
        par_change = [0.9 0.8 1.1 1.2];
    % CTEPH treatment - surgical intervention
    case 5
        par_change = [1.0 0.2 1.0 1.8];
    case 6
        par_change = [1.0 0.4 1.0 1.6];
    case 7
        par_change = [1.0 0.6 1.0 1.4];
    % CTEPH treatment - surgery & vasodilator
    case 8
        par_change = [0.8 0.15 1.2 1.85];
    case 9
        par_change = [0.8 0.35 1.2 1.65];
    case 10
        par_change = [0.8 0.55 1.2 1.45];
end

pars(INDMAP) = pars(INDMAP).*par_change';
%%


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
Trrv = Tcrv + pars(25);
% right ventricle relaxation

Tclv = Tcrv.*(0.95); % left ventricle contraction MITCHEL SCALE!
Trlv = Trrv; % left ventricle relaxation
    
% Initialize vectors for pressure solutions
psaS  = [];  % systemic arteries
psvS  = [];  % systemic veins
ppaS  = [];  % pulmonary arteries
ppvS  = [];  % pulmonary veins
praS  = [];  % right atrium
plaS  = [];  % left atrium
prvS  = [];  % right ventricle
plvS  = [];  % left ventricle

ElaS  = [];  % Elastance left atrium

% Initialize vectors for volume solutions
VsaS  = [];  % systemic arteries
VsvS  = [];  % systemic veins
VpaS  = [];  % pulmonary arteries
VpvS  = [];  % pulmonary veins
VraS  = [];  % right atrium
VlaS  = [];  % left atrium
VrvS  = [];  % right ventricle
VlvS  = [];  % left ventricle

% Initialize vectors for flow solutions
qsS   = []; % systemic periphery
qpS   = []; % pulmonary periphery
qmvaS = []; % flow between left atrium  and ventricle (mitral valve)
qavaS = []; % flow out of aortic valve
qtvaS = []; % flow between right atrium and ventricle (tricuspid valve)
qpvaS = []; % flow out of pulmonary arterial valve
qsvS  = []; % flow into right atrium (systemic veins)
qpvS  = []; % flow into left atrium (pulmonary veins)
    
% Initialize CO vectors
COsS  = [];
COpS  = [];
tCOS  = [];

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
   qava = zeros(size(qs));
   qpva = zeros(size(qs));
   qtva = zeros(size(qs));
   qsv  = zeros(size(qs));
   qpv  = zeros(size(qs));
   for j = 1:length(tdc)
      if pla(j) > plv(j) % Mitral valve
         qmva(j) = (pla(j)-plv(j))/Rmva; % valve open
      end
    
      if plv(j) > psa(j) % Aortic valve
         qava(j) = (plv(j)-psa(j))/Rava; % valve open
      end
    
      if prv(j) > ppa(j) % flow through pulmonary valve
        qpva(j) = (prv(j)-ppa(j))/Rpva;
      end
    
      if pra(j) > prv(j)  % Tricuspid valve
         qtva(j) = (pra(j)-prv(j))/Rtva; % valve open
      end

      if psv(j) > pra(j)
        qsv(j) = (psv(j)-pra(j))/Rsv;
      end

      if ppv(j) > pla(j)
       qpv(j) = (ppv(j)-pla(j))/Rpv;
      end
   end;
    
    % save solution
    psaS  = [psaS  psa(1:end-1)'];  % systemic arteries
    psvS  = [psvS  psv(1:end-1)'];  % systemic veins
    ppaS  = [ppaS  ppa(1:end-1)'];  % pulmonary arteries
    ppvS  = [ppvS  ppv(1:end-1)'];  % puylmonary veins
    praS  = [praS  pra(1:end-1)'];  % right atrium
    plaS  = [plaS  pla(1:end-1)'];  % left atrium
    prvS  = [prvS  prv(1:end-1)'];  % right ventricle
    plvS  = [plvS  plv(1:end-1)'];  % left ventricle

    ElaS  = [ElaS  Ela(1:end-1)];

    VsaS   = [VsaS  Vsa(1:end-1)']; % systemic arteries
    VsvS   = [VsvS  Vsv(1:end-1)']; % systemic veins
    VpaS   = [VpaS  Vpa(1:end-1)']; % pulmonary arteries
    VpvS   = [VpvS  Vpv(1:end-1)']; % puylmonary veins
    VraS   = [VraS  Vra(1:end-1)']; % right atrium
    VlaS   = [VlaS  Vla(1:end-1)']; % left atrium
    VrvS   = [VrvS  Vrv(1:end-1)']; % right ventricle
    VlvS   = [VlvS  Vlv(1:end-1)']; % left ventricle

    qpS    = [qpS   qp(1:end-1)'];   % systemic periphery
    qsS    = [qsS   qs(1:end-1)'];   % pulmonary periphery
    qavaS  = [qavaS qava(1:end-1)'];
    qpvaS  = [qpvaS qpva(1:end-1)'];
    qmvaS  = [qmvaS qmva(1:end-1)']; % flow between left atrium  and ventricle (mitral valve)
    qtvaS  = [qtvaS qtva(1:end-1)']; % flow between right atrium and ventricle (tricuspid valve)
    qsvS   = [qsvS  qsv(1:end-1)'];  % flow into right atrium (systemic veins)
    qpvS   = [qpvS  qpv(1:end-1)'];  % flow into left atrium (pulmonary veins)

    COsS   = [COsS COs];
    COpS   = [COpS COp];
    tCOS   = [tCOS tdc(1)];        
    Init = [Vsa(end) Vsv(end) Vpa(end) Vpv(end) Vra(end) Vla(end) Vrv(end) Vlv(end)];% qava(end) qpva(end)];
    if i < NC
       k1 = k2; %sets last index of this loop as first index for next loop
       k2 = round(T/dt)+k1;
    end
end

 % Save last time point
psaS  = [psaS  psa(end)];  % systemic arteries
psvS  = [psvS  psv(end)];  % systemic veins
ppaS  = [ppaS  ppa(end)];  % pulmonary arteries
ppvS  = [ppvS  ppv(end)];  % puylmonary veins
praS  = [praS  pra(end)];  % right atrium
plaS  = [plaS  pla(end)];  % left atrium
prvS  = [prvS  prv(end)];  % right ventricle
plvS  = [plvS  plv(end)];  % left ventricle

ElaS  = [ElaS  Ela(end)];

VsaS   = [VsaS  Vsa(end)]; % systemic arteries
VsvS   = [VsvS  Vsv(end)]; % systemic veins
VpaS   = [VpaS  Vpa(end)]; % pulmonary arteries
VpvS   = [VpvS  Vpv(end)]; % puylmonary veins
VraS   = [VraS  Vra(end)]; % right atrium
VlaS   = [VlaS  Vla(end)]; % left atrium
VrvS   = [VrvS  Vrv(end)]; % right ventricle
VlvS   = [VlvS  Vlv(end)]; % left ventricle

qpS   = [qpS   qp(end)];   % systemic periphery
qsS   = [qsS   qs(end)];   % pulmonary periphery
qavaS = [qavaS qava(end)];
qpvaS = [qpvaS qpva(end)];
qmvaS = [qmvaS qmva(end)]; % flow between left atrium  and ventricle (mitral valve)
qtvaS = [qtvaS qtva(end)]; % flow between right atrium and ventricle (tricuspid valve)
qsvS  = [qsvS  qsv(end)];  % flow into right atrium (systemic veins)
qpvS  = [qpvS  qpv(end)];  % flow into left atrium (pulmonary veins)



Xlimits =  [tdc(1) tdc(end)];

figure(2); clf;
subplot(3,3,1); hold on;
plot(tdc,pla,'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('pla (mmHg)');
xlim(Xlimits);
hold off;

subplot(3,3,2); hold on;
plot(tdc,plv,'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('plv (mmHg)');
xlim(Xlimits);
hold off;

subplot(3,3,3); hold on;
plot(tdc,psa,'r','LineWidth',3.2);
set(gca,'FontSize',24);
plot(tdc,data.pSAM*ones(size(tdc)),'--k','LineWidth',3.2);
plot(tdc,data.pSAm*ones(size(tdc)),'--k','LineWidth',3.2);
ylabel('psa (mmHg)');
xlim(Xlimits);
hold off;

subplot(3,3,4); hold on;
plot(tdc,psv,'r','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('psv (mmHg)');
xlim(Xlimits);
hold off;

subplot(3,3,5); hold on;
plot(tdc,pra,'r','LineWidth',3.2);
plot(tdc, max(data.pRA)*ones(size(tdc)), '--k', 'LineWidth', 2);
plot(tdc, min(data.pRA)*ones(size(tdc)), '--k', 'LineWidth', 2);
plot(tdc,data.pRA,'k','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('pra (mmHg)');
xlim(Xlimits);
hold off;

subplot(3,3,6); hold on;
plot(tdc,prv,'r','LineWidth',3.2);
plot(tdc, max(data.pRV)*ones(size(tdc)), '--k', 'LineWidth', 2);
plot(tdc, min(data.pRV)*ones(size(tdc)), '--k', 'LineWidth', 2);
plot(tdc,data.pRV,'k','LineWidth',3.2);
xlim(Xlimits);
set(gca,'FontSize',24);
ylabel('prv (mmHg)');
hold off;


subplot(3,3,7);hold on;
plot(tdc,ppa,'r','LineWidth',3.2);
plot(tdc, max(data.pPA)*ones(size(tdc)), '--k', 'LineWidth', 2);
plot(tdc, min(data.pPA)*ones(size(tdc)), '--k', 'LineWidth', 2);
plot(tdc,data.pPA,'k','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('ppa (mmHg)');
xlabel('time (s)');
xlim(Xlimits);
hold off;
    
    
subplot(3,3,8);hold on;
plot(tdc,mean(ppv)*ones(size(tdc)),'r','LineWidth',3.2);
plot(tdc, min(data.pPW)*ones(size(tdc)), '--k', 'LineWidth', 2);
set(gca,'FontSize',24);
ylabel('ppv (mmHg)');
xlabel('time (s)')
ylim([mean(ppv)*.9 , mean(ppv)*1.1])
xlim(Xlimits);
hold off;

subplot(3,3,9);hold on;
plot(tdc,ones(size(tdc)).*COsS(end),'r','LineWidth',3.2);
plot(tdc,ones(size(tdc)).*data.CO,'--k','LineWidth',3.2);
set(gca,'FontSize',24);
ylabel('CO (mL/sec)');
xlabel('time (s)')
ylim([round(data.CO)-2 round(data.CO)+2])
hold off;

h = figure(2);
h.Position = [379   778   951   552]; legend off
end 


