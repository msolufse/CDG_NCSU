%%
% Function that will take a geometry in and change everything to have a
% stenosis in one of the vessels.

function [conn] = make_stenosis_2(conn,sten_L)
% Take connectivity file and add blockage
% DOES NOT WORK FOR TRIFURCATIONS AT THIS TIME
[ids,~,~] = find(conn(:,1)==sten_L);
conn(ids,4) = -1;
end