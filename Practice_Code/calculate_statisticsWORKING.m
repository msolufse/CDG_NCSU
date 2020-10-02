function stats = calculate_statisticsWORKING(r,c,data)
% r - row, c - column/subset, 

load('PostOptP4_WSubsets.mat')
history = PostHistory;
[~,Init,~,~,~] = load_global_PH(data);
data.Init = Init;

%% These are examples of structs with multiple cells of information within them %%
Pars = IND{r,c};
ppa  = history{r,c}.ppa;
prv  = history{r,c}.prv;
pra  = history{r,c}.pra;

%%
model = [ppa prv pra]; 


dPrv = data.pRV;
dPpa = data.pPA;
dPra = data.pRA;
DATA = [dPpa dPrv dPra];

Np = length(Pars); % Number of parameters
% Nd = length(data(1,:)); % Number of data points

% R2 calulcation
for i = 1:3
    ybar(i)  = mean(DATA(:,i));
    SStot(i) = sum((DATA(:,i)-ybar(i)).^2);      % Variance of the data
    SSmod(i) = sum((DATA(:,i)-model(:,i)).^2);   % Unexplained error (model mismatch)
    R2(i)    = 1 - SSmod(i)./SStot(i);           % R2 value
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compute R^2 value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
mdl   = fitlm(model(:,1),DATA(:,1));
R22(1) = mdl.Rsquared.Ordinary;
 
mdl   = fitlm(model(:,2),DATA(:,2));
R22(2) = mdl.Rsquared.Ordinary;
 
mdl   = fitlm(model(:,3),DATA(:,3));
R22(3) = mdl.Rsquared.Ordinary;


%SSmod for AIC is rout from CVmodel
% call CV model with a specific residual (eg 8)
pars = history{r,c}.pars;
data.res_flag = 4;
[rout,J,Lts]= CVmodel(pars,data);
Nd = length(rout); % Number of data points

% Res = ones(length(rout),1);
% Res(1:Lts(3)) = rout(1:Lts(3))*sqrt(Lts(3));
% Res(Lts(3)+1:end) = rout(Lts(3)+1:end)*sqrt(Lts(4));
Res = rout*sqrt((Lts(1)+Lts(2)));
J2 = Res'*Res; % should be a scalar


% sig2 = J2./(Nd - Np); %variance
% likelihood = exp(-J2/(2*sig2))./(2*pi*sig2)^(Nd/2); %Likelihood


% Compute AIC
AIC = 2*Np + 2*log(J); %added 2* to the log likelihood  %changed - to + 

% Compute AICc: corrected for small sample size
AICc = AIC + 2*(Np^2 + Np)./(Nd-Np-1);

% BIC: Bayesian information criteria; penalizes differently for # params
% Note that BIC works better when comparing multiple models if we can get
% very close to the measured data

BIC = 2*log(Nd)*Np + 2*log(J); % changed - to + 


%% These are more structs %%

stats.AIC  = AIC;
stats.AICc = AICc;
stats.BIC  = BIC;
stats.R2   = R2;
stats.J    = J;
stats.R22  = R22;

end


