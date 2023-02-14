% This function initializes the parameters for the model and sets initial
% values for the variables.
function [pars,Vunstr,Init,data] = load_global_CTR
data.INDMAP   = 1:25; % Parameter index (used in optimization)
data.ODE_TOL  = 1e-8; % Absolute ODE tolerance
data.REL_TOL  = 1e-8; % Relative ODE tolerance
data.DIFF_INC = 1e-4; % Difference increment
data.Height   = 180;  % Average female height
data.Weight   = 75;   % Average female weight
data.Gender   = 2;    % Female

%% Total blood volume (mL)
BSA = sqrt(data.Height*data.Weight/3600);
data.BSA    = BSA; 

if data.Gender == 2
  TotalVol = (3.47*BSA - 1.954)*1000; % Female
elseif data.Gender == 1
  TotalVol = (3.29*BSA - 1.229)*1000; % Male
end

%% Cardiac Output (mL/s)
data.CO = 5;            % Reference cardiac ouput
qtot = data.CO*1000/60; % Total flow (ml/s)
data.T = 1;             % Cardiac Cycle Length
data.NC = 30;           % Number of heart beats to solve the model
data.dt = 0.001;        % Time stepping for model

% Construct a time vector
data.td = [0:data.dt:data.T*data.NC];

%% Stroke Volume
Vstroke   = qtot*data.T; % (mL/s * s/beat) = mL/beat)
data.Vstr = Vstroke;     

%% Blood Pressure

% Systemic arteries
psa    = 93;       % mean systemic artery pressure
pSsa   = 120;      % systolic systemic artery pressure
pDsa   = 80;       % systolic systemic artery pressure

% Pulmonary arteries
ppa    = 12;       % mean pulmonary artery pressure  
pSpa   = 21;       % systolic pulmonary artery pressure
pDpa   = 8;        % systolic pulmonary artery pressure

% Veins
ppv    = 5;        % mean pulmonary venous pressure
psv    = 8;        % mean systemic venous pressure

%% The heart
% Atria
pmra   = 2;       % minimum right atrial pressure
pMra   = 12;      % maximum right atrial pressure
pmla   = 2;       % minimum left atrial pressure
pMla   = 13;      % maximum left atrial pressure             

% Ventricles
pMrv   = 1.01*pSpa;   % maximum right ventricular pressure  
pmrv   = 3;           % minimum right ventricular pressure 
pMlv   = 1.01*pSsa;   % maximum left ventricular pressure         
pmlv   = 8;           % minimum left ventricular pressure

% Blood Volumes
VraM  = BSA*30;      %REF
Vram  = BSA*13;      %REF

VlaM  = BSA*30;      %REF
Vlam  = BSA*13;      %REF

VraMF = VraM/TotalVol; %REF
VlaMF = VlaM/TotalVol; %REF

VrvM = 78*BSA;         %REF
VlvM = 78*BSA;         %REF
Vlvm = VlvM-Vstroke;   %REF    %left ventricular min volume
Vrvm = VrvM-Vstroke;   %REF    %right ventricular max

VrvMF = VrvM/TotalVol; %REF 
VlvMF = VlvM/TotalVol; %REF

HVF  = VraMF+VlaMF+VrvMF+VlvMF;

% Volume distribution 
Vsa = 0.13*TotalVol;       % REF Arterial volume (13%)
Vsv = (0.72-HVF)*TotalVol; % REF Venous volume (~65%)  
Vpa = 0.03*TotalVol;       % REF Pulmonary artery (3%)
Vpv = 0.11*TotalVol;       % REF Pulmonary vein volume (11%)
 
% Stressed Volumes (from Beneken)
VsaU = Vsa*(1-0.27);  % Stressed volume 27% of arterial volume
VsvU = Vsv*(1-0.08);  % Stressed volume 7.5% of venous volume
VpaU = Vpa*(1-0.58);  % Stressed volume pulmonary artery 
VpvU = Vpv*(1-0.11);  % Streseed volume pulmonary vein 
VaU  = 10;           
VvU  = 5;
Vunstr = [VsaU VsvU VpaU VpvU VaU VvU];

data.V = Vunstr;

data.VlvM = VlvM;
data.VrvM = VrvM;
data.VlaM = VlaM;
data.VraM = VraM;
data.Vlam = Vlam;
data.Vram = Vram;

data.RaVstr = VraM-Vram;
data.LaVstr = VlaM-Vlam;

%% Resistences
% Peripheral resistances (Ohm's law)
Rs   = (pDsa-psv)/qtot;  % systemic resistance 
Rp   = (pDpa-ppv)/qtot;  % pulmonary resistance

% Valve resistances
Rava = (pMlv-pSsa)/qtot; % resistance aortic valve
Rpva = (pMrv-pSpa)/qtot; % resistance pulmonary valve 

Rmva = 0.25./qtot;  %REF
Rtva = 0.25./qtot;  %REF

% Venous resistances
Rpv = (ppv-pmla)/qtot;       % resistance of pulmonary vein
Rsv = (psv-pmra)/qtot;       % resistance of vena cava

%% Compliances
Csa = (Vsa-VsaU)/psa;       % systemic artery compliance
Csv = (Vsv-VsvU)/psv;       % systemic venous compliance  
Cpa = (Vpa-VpaU)/ppa;       % pulmonary artery compliance
Cpv = (Vpv-VpvU)/ppv;       % pulmonary vein compliance

%% Heart parameters
% Elastance parameters atria 
Emra   = pmra/(VraM-VaU); %min right atrium 
EMra   = pMra/(Vram-VaU); %max right atrium 
Emla   = pmla/(VlaM-VaU); %min left atrium
EMla   = pMla/(Vlam-VaU); %max left atrium

% Elastance parameters ventricles
Emrv   = pmrv/(VrvM-VvU); %min right ventricle
EMrv   = pMrv/(Vrvm-VvU); %max right ventricle    
Emlv   = pmlv/(VlvM-VvU); %min left ventricle
EMlv   = pMlv/(Vlvm-VvU); %max left ventricle

% Timing heart contraction 
% Aria
Trra = 0.10;  %REF 
tcra = 0.80;  %REF 
Tcra = 0.97;  %REF

% Ventricle
Tcrv = 0.3;  %REF
Trrv = 0.51; %REF

%% Initial conditions
Init = [Vsa Vsv Vpa Vpv VraM VlaM VrvM VlvM];
data.Init = Init;
%% Parameter vector
% LIST ALL PARAMETERS I A TABLE 
x0 = [Rs Rp Rava Rmva Rpva Rtva Rpv Rsv ...  % 1-8
      Csa Csv Cpa Cpv ...                    % 9-12
      EMra Emra EMla Emla ...                % 13-16
      EMrv Emrv EMlv Emlv ...                % 17-20
      Trra tcra Tcra ...                     % 21-23
      Tcrv Trrv]';                           % 24-25
pars = log(x0);

