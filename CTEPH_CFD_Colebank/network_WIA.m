% Get peak WIA values throughout the network

function out = network_WIA(t,p,q,a,r,stiff,plot_ID)
% Define constants
rho = 1.055;
n = length(p);
p = p.*1333.22;
u = q./a;
dp = diff(p);
du = diff(u);
dt = t(2)-t(1);

% NOTE: W/m^2 is 0.001 g/s^3
Wconv = 0.001;

% Fudge factor: C was not nondimensionalized correctly
c = wave_speed(a,r,stiff);

rho_cW = rho.*c(1:end-1,:);

dp_fwd = 0.5.*(dp+rho_cW.*du);
dp_bwd = 0.5.*(dp-rho_cW.*du);

P_fwd = fwd_bwd(dp_fwd,p);
P_bwd = fwd_bwd(dp_bwd,p);

du_fwd = 0.5.*(du+(dp./rho_cW));
du_bwd = 0.5.*(du-(dp./rho_cW));

U_fwd = fwd_bwd(du_fwd,u);
U_bwd = fwd_bwd(du_bwd,u);

WI_fwd = Wconv.*(dp_fwd.*du_fwd)./dt^2;
WI_bwd = Wconv.*(dp_bwd.*du_bwd)./dt^2;

figure(plot_ID);clf; 
semilogy(max(WI_fwd),'-ko');
hold on;
semilogy(abs(min(WI_bwd)),'-ro');

end

function c = wave_speed(a,r,fr)
a0 = pi.*r.^2;
a0 = a0';
fR = fr(r)';
c =  sqrt(0.5.*fR.*sqrt(a./a0));
% c =  0.5*fr(r)*sqrt(a/a0); % use if fr is a function
end

function out = fwd_bwd(dX,X)
n = size(dX,1);
out = zeros(n,size(dX,2));
out(1,:) = dX(1,:);
for i=2:n
    out(i,:) = out(i-1,:)+dX(i,:);
end
out = out + X(1,:);
end

