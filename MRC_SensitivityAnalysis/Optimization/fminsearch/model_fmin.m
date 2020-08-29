function [ss,sol,rout] = model_fmin(pars)

global data

xdata = data.xdata;
ydata = data.ydata;

sol = model_wrap(pars);

rout = (sol(:,2) - ydata);
ss = rout'*rout;
