%--------------------------------------------------------------------------
%Define ODEs, using parameters from load_global.m
%--------------------------------------------------------------------------

function ydot = modelBasic(t,y,pars)

y1  = y(1);
y2  = y(2);

gamma = pars(1);	k = pars(2);	r = pars(3);	delta = pars(4);

dy1 = delta*1000 - delta*y1 - gamma*k*y2*y1;
dy2 = gamma*k*y2*y1 - (r + delta)*y2;

ydot = [dy1; dy2];