function [Ehr0, R1, R2, CT, L] = ND_Par (E, h, r0, Qmean, Pmean, alpha, tau, Inductance)

% This function comutes and nondimensionalise the arterial stiffenss and
% winddkessel parameters.

rho = 1.055;
g = 981;
Lr = 1.0;
qc = 10*Lr^2;
conv = 1332.22;

tc = tau*qc/Lr^3;

CO = Qmean/qc;
mP = Pmean*conv/rho/g/Lr;

RT = mP/CO;
R1 = alpha*RT;
R2 = RT - R1;
CT = tc/RT;

L = Inductance*conv*qc^2/rho/g/Lr^4;

Ehr0 = conv*E*h/r0;

end