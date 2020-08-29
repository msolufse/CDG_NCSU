function DriverBasic_intervals

global ALLPARS INDMAP Init

load optResults.mat %parsNM, parsLM, histout, jachist, rout, and sc
%load optResultsWeighted.mat %r = (ymodel-ydata)/mean(ydata)

parsopt = parsLM;   
ALLPARS = parsopt; 


[x0,Init,~,~] = load_global;

load SIR_data.mat
ydata = obs;
xdata = (0:0.1:6)';

INDMAP = 2:4; %FIND A WAY TO GET THIS FROM COVARIANCE.M AUTOMATICALLY!

% Integrate the system using ode45.m
sol = model_sol(parsopt,xdata,Init);
Isol = sol(2,:)';

figure(1); hold on
h = plot(xdata,Isol,xdata,ydata,'*');
set(h,'linewidth',4);
set(h,'Markersize',4);
set(gca,'Fontsize',24);
xlabel('Time')
ylabel('Number of infections')
grid on

Jac = jachist{end};
Jac = bsxfun(@rdivide, Jac,diag(sc)'); %not really necessary since bounds 
                                       %are 0 and 1 here, but need otherwise

%-----------------------Parameter confidence intervals---------------------
[N,M] = size(Jac);
sig = 1/(N - M) * (rout'*rout); %mse
COV = sig*(inv(Jac'*Jac)); %

CI = nlparci(parsopt(INDMAP),rout,'covar',COV)
parsopt(INDMAP) - CI(:,1)
pause

%------------------------------Model intervals-----------------------------
t   = linspace(0,6,200)';

w = ones(1,length(t))';
%w = ones(1,length(t))'/((mean(obs))^2); %use if 
                                        %rout = (ymodel - ydata)/mean(ydata)

%newlsq2
[ypredC1,deltaC1] = nlpredci(@model_wrap,t,parsopt(INDMAP),rout,'Covar',COV,'MSE',sig,'Weights',w);
[ypredP1,deltaP1] = nlpredci(@model_wrap,t,parsopt(INDMAP),rout,'Covar',COV,'MSE',sig,'Weights',w,'PredOpt','observation');

%Exactly the same for the two optimizers, so only plot one!
figure(1); hold on
m = plot(t,ypredC1 - deltaC1,'b--',t,ypredC1 + deltaC1,'b--');
set(m,'linewidth',2);
n = plot(t,ypredP1 - deltaP1,'g--',t,ypredP1 + deltaP1,'g--');
set(n,'linewidth',2);
set(gca,'Fontsize',24);
legend([m(1),n(1)],'95% Confidence Interval','95% Prediction Interval')