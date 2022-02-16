function [p,q,a,c] = STmodel_V2(pars,sc)
if ~isempty(sc)
    pars = [pars(1:end-2).*sc pars(end-1:end)];
end
% Use rmin as a parameter
tic
flag = unix(sprintf('./sor06 %f %f %f %f %f %f %f %f %f %f %f %d %d',pars));
toc

data = load('pu_ALL.2d');
if flag ~= 1 && ~isempty(data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Run this to get all the vessels
    if ~isempty(data)
        [~,~,p,q,a,c] = gnuplot(data);
    end
else
    flag = unix(sprintf('./sor06 %f %f %f %f %f %f %f %f %f %f %f %d %d',pars));
    data = load('pu_ALL.2d');
    [~,~,p,q,a,c] = gnuplot(data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
