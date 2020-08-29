%--------------------------------------------------------------------------
%Sum of squares function
%Uses model_wrap.m to fix some parameters and let the others vary (INDMAP)
%--------------------------------------------------------------------------
function [ss,sol,rout] = model_res(pars,data)

xdata  = data.xdata;
ydata  = data.ydata;
Init   = data.Init;

sol    = model_wrap(xdata,pars,Init);

rout = (sol(:,2) - ydata);
% rout = (sol(:,2) - ydata)/mean(ydata);
ss = rout'*rout;
 