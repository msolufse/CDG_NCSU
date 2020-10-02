function [Q] = InFlow(Qmin,Qmax, t0,Tm,tau,T )

% This function creates an inflow profile to be imposed at the vessel inlet

N = 8192;
 
deltaT = (T-t0)/N;

t = t0:deltaT:T;


Q = zeros(1,N+1);

for i = 1:N+1
    
    ti = i*deltaT;
    
    if (t0<=ti)&&(ti<=Tm)
        Q(i) = 0.5*(Qmax - Qmin)*(1 - cos(pi*ti/Tm)) + Qmin;
    elseif (Tm<=ti)&&(ti<=tau+Tm)
        Q(i) = 0.5*(Qmax - Qmin)*(1 + cos(pi*(ti-Tm)/tau)) + Qmin;
    elseif (tau+Tm<=ti)&&(ti<=T)
        Q(i) = Qmin;

    end
    
end

fname = strcat('Qin_',num2str(N),'.dat');

dlmwrite(fname,Q');

figure(1);
plot(t,Q,'linewidth',2);
set(gca, 'fontsize',20,'xlim',[t0 T]);grid on;
xlabel 'Time (s)';
ylabel 'Inflow (ml/s)';
title 'Inflow boundary condition';


end

