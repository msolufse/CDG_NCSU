% This function initializes the parameters for the model and sets initial
% values for the variables.
function [pars, Init,low,hi,qtot] = load_global_PH(data)


global DIFF_INC ODE_TOL REL_TOL 

ODE_TOL  = 1e-8;
REL_TOL  = 1e-8;
DIFF_INC = sqrt(ODE_TOL);
%DIFF_INC = 1e-2;

%% Patient-specific BSA and total blood volume
BSA = sqrt(data.Height*data.Weight/3600); %Body surface area (Shoemaker W (1989) Fluids)
if     data.Gender == 2
  TotalVol = (3.47*BSA - 1.954)*1000; % Female
elseif data.Gender == 1
  TotalVol = (3.29*BSA - 1.229)*1000; % Male
end

%% Patient-specific Cardiac Output
% qtot = data.CO*1000/60;  % data from thermodilution
qtot = TotalVol/60;        % predicting volume from BSA

%% Patient-specific Pressures 
% Systemic arteries
psa    = data.pSAa;                   % mean systemic artery pressure
pSsa   = data.pSAM;                   % systolic systemic artery pressure
pDsa   = data.pSAm;                   % diastolic systemic artery pressure 

% Pulmonary arteries
ppa    = mean(data.pPA);                % mean pulmonary artery pressure (area under curve) 
pSpa   = max(data.pPA);                 % systolic pulmonary artery pressure
pDpa   = data.pPA(1);                   % diastolic pulmonary artery pressure

% Veins
ppv    = min(data.pPW);
psv    = 1.01*min(data.pRA);            % systemic veins 

% Atria
pmra   = min(data.pRA);
pMra   = max(data.pRA);
pmla   = 0.95*ppv;                       % min left atrial pressure
pMla   = pmla + 5; %Was 10, too big MJC  % Find literature value

% Ventricles
pMrv   = max(data.pRV);                 % max right ventricular pressure
pmrv   = min(data.pRV);                 % min right ventricular pressure
pMlv   = 1.01*pSsa;                     % max left ventricular pressure
pmlv   = 0.97*pMla;                     % min left ventricular pressure

%% Patient-specific Volumes
% Volume distribution
VTsa = 0.13*TotalVol; % Arterial volume (13%)
VTsv = 0.65*TotalVol; % Venous volume (65%)
VTpa = 0.03*TotalVol; % Pulmonary artery (3%)
VTpv = 0.11*TotalVol; % Pulmonary vein volume (11%)
 
% Stressed Volumes (from Beneken)
VSsa = VTsa*0.27;  % Stressed volume 27% of arterial volume
VSsv = VTsv*0.075; % Stressed volume 7.5% of venous volume
VSpa = VTpa*0.60;  % Stressed volume pulmonary artery was (RECALC)
VSpv = VTpv*0.11;  % Streseed volume pulmonary vein (RECALC)

% Ejection fraction 
EF_V = 0.60; % Ventricular ejection fraction
EF_A = 0.47; % Atrial ejection fraction

% Atrial volume
Vda  = 5;
VMra = 0.015*TotalVol-Vda;  %right atrial max minus dead space volume
Vmra = VMra*(1-EF_A);       %right atrial min
VMla = VMra;                %left atrial max volume (stressed)
Vmla = Vmra;                %right atrial min

% Ventricular volume
Vdv  = 10;
VMrv = 0.025*TotalVol-Vdv;  %right ventricular max (stressed)
Vmrv = VMrv*(1-EF_V);       %right ventricular min
VMlv = VMrv;                %left ventricular volume (stressed)
Vmlv = Vmrv;                %left ventricular min

%% Resistences
% Peripheral resistances (Ohm's law)
Rs  = (psa-psv)/qtot;   % systemic resistance
Rp  = (ppa-ppv)/qtot;   % pulmonary resistance

% Valve resistances
Rava = (pMlv-pSsa)/qtot; % resistance aortic valve
Rpva = (pMrv-pSpa)/qtot; % resistance pulmonary valve (Data too big a jump - not simultaneously measured)

Rmva = (pMla-pmlv)/qtot; % resistance mitral valve
Rtva = 0.055; %(pMra-pmrv)/qtot  % resistance tricuspid valve (Data too big a jump - not simultaneously measured)

% Venous resistances
Rpv = (ppv-pmla)/qtot;    % resistance of pulmonary vein
Rsv = (psv-pmra)/qtot;    % resistance of vena cava

%% Compliances
Csa = VSsa/pSsa;   % systemic artery compliance
Csv = VSsv/psv;    % systemic venous compliance  
Cpa = VSpa/pSpa;   % pulmonary artery compliance
Cpv = VSpv/ppv;    % pulmonary vein compliance

%% Heart parameters
% Elastance parameters atria
Emra   = pmra/VMra; %min right atrium
EMra   = pMra/Vmra; %max right atrium 
Emla   = pmla/VMla; %min left atrium
EMla   = pMla/Vmla; %max left atrium

% Elastance parameters ventricles
Emrv   = pmrv/VMrv; %min right ventricle
EMrv   = pMrv/Vmrv; %max right ventricle    
Emlv   = pmlv/VMlv; %min left ventricle
EMlv   = pMlv/Vmlv; %max left ventricle

tcra = data.tcra;
Trra = data.Trra;
Tcra = data.Tcra;
Tcrv = data.Tcrv;
Trrv = data.Trrv;

%% Initial conditions
Init = [VSsa VSsv VSpa VSpv Vmra Vmla VMrv VMlv];% Steady state

%% Parameter vector
x0 = [Rs Rp Rava Rmva Rpva Rtva Rpv Rsv ...  % 1-8
      Csa Csv Cpa Cpv ...                    % 9-12
      EMra Emra EMla Emla ...                % 13-16
      EMrv Emrv EMlv Emlv ...                % 17-20
      Trra tcra Tcra ...                     % 21-23
      Tcrv Trrv]';                           % 24-25

pars = log(x0);

hi  = pars + log(4);
low = pars - log(10);

hi(21:25)  = pars(21:25) + log(1.5);
low(21:25) = pars(21:25) - log(1.5);

