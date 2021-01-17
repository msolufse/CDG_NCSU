% This function initializes the parameters for the model and sets initial
% values for the variables.
function [pars, Init,hi,low] = load_global_CONTROL %call number

global DIFF_INC ODE_TOL 

ODE_TOL  = 1e-8;
DIFF_INC = 1e-4;

H  = 164;
W  = 70;
G  = 1;
CO = 6;
% Define subject specific quantities
BSA = sqrt(H*W/3600);% load in   %Body surface area

if G == 2
  TotalVol = (3.47*BSA - 1.954)*1000; % Female
elseif G == 1
  TotalVol = (3.29*BSA - 1.229)*1000; % Male
end
% qtot = CO*1000/60; %TotalVol/60 67 Assume blood circulates in one min
qtot = TotalVol/60;

% Flows (related to subject)
qs = qtot; % Upper body flow arteries --> veins
qp   = qtot; % Upper body arteries --> lower body arteries

% Pressures (related to subject) 
% MJC: All pressure values are obtained using a healthy, textbook
% definition. These values can all be found in the textbook
% "Medical Physiology" by Boron and Boulpaep
psa    = 95; % mean systemic artery pressure
pSsa   = 120;%120; % systolic systemic artery pressure               % 
pDsa   = 80;%80; % diastolic systemic artery pressure 
psv    = 15;  % mean systemic vein pressure
ppa    = 15; % mean pulmonary artery pressure 
pSpa   = 15;%20;  % systolic pulmonary artery pressure
pDpa   = 8;  % diastolic pulmonary artery pressure
ppv    = 10;  % mean pulmonary vein pressure
pMla   = 12; % maximum left atrial pressure           
pmla   = 8;  % minimum left atrial pressure
pMlv   = 130;%125;% maximum left ventricular pressure
pmlv   = 10; % minimum left ventricular pressure
pMra   = 8;  % maximum right atrial pressure 
pmra   = 3;  % minimum right atrial pressure
pMrv   = 18;%25; % maximum right ventricular pressure
pmrv   = 4;  % minimum right ventricular pressure

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
Rs  = (pSsa-psv)/qtot;   % systemic resistance
Rp  = (pSpa-ppv)/qtot;   % pulmonary resistance

% Valve resistances
Rava = (pMlv-pSsa)/qtot; % resistance aortic valve
Rpva = (pMrv-pSpa)/qtot; % resistance pulmonary valve (Data too big a jump - not simultaneously measured)

Rmva = (pMla-pmlv)/qtot; % resistance mitral valve
Rtva = (pMra-pmrv)/qtot;  % resistance tricuspid valve (Data too big a jump - not simultaneously measured)

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

Trra = 0.25;
tcra = 0.40;
Tcra = 0.20;
Tcrv = 0.30;
Trrv = 0.20;


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

% Optimizing timing parameters check upper and lower bounds
% hi(21:25)  = pars + log(2);
% low(21:25) = pars - log(2);


  
  