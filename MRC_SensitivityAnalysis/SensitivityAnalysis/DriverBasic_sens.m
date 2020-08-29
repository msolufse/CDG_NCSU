%--------------------------------------------------------------------------
%Computes the sensitivity matrix dy/dpars.
%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set 
%pars = log(pars) in DriverBasic_sens
%Plots ranked sensitivities 
%--------------------------------------------------------------------------

function DriverBasic_sens

global ydata xdata

load SIR_data.mat %obs

[pars,Init] = load_global;
xdata = (0:0.01:6)';
ydata = obs;

%senseq finds the non-weighted sensitivities
sens = senseq(pars,xdata,Init);

% ranked classical sensitivities
[M,N] = size(sens);
for i = 1:N
  sens_norm(i)=norm(sens(:,i),2);
end

[Rsens,Isens] = sort(sens_norm,'descend');
display([Isens]);

%Ranked sensitivities
figure(1);clf; 
h=semilogy(Rsens./max(Rsens),'x');
set(h,'linewidth',4);
set(h,'Markersize',24);
set(gca,'Fontsize',24);
grid on;
ylim([0.025 1]);
ylabel('Sensitivities');
xlabel('Parameters')
print -depsc2 RankedSensitivities.eps

%Plots time-varying sensitivities
figure(3); clf;
h=plot(xdata,sens(:,1), ...
       xdata,sens(:,2), ...
       xdata,sens(:,3), ...
       xdata,sens(:,4));
set(h,'Linewidth',4);
set(gca,'Fontsize',24);
xlabel('Time');
ylabel('Sensitivities');
grid on;
print -depsc2 TimevaryingSensitivities.eps

save Sens.mat sens Rsens Isens;