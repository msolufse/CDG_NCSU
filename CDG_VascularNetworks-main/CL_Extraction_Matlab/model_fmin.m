function [ss,sol,rout] = model_fmin(pars)

global datac

xdc = datac.xdc;
rdc = datac.rdc;

k1 = pars(1);
k2 = pars(2);
k3 = pars(3);

rad = k1*exp(-k2*xdc)+k3;

rout = (rad - rdc)';
ss = rout'*rout;

