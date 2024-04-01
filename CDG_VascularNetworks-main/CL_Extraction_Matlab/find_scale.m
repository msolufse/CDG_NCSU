function [ alpha ] = find_scale(scaling, vessel)
%% Takes parent vessel, scales it, and applies scaling factor to all other vessels
z_temp = vessel(:,3);
if z_temp(1) > z_temp(end)
    vessel = flipud(vessel);
end
x = vessel(:,1);
y = vessel(:,2);
z = vessel(:,3);
r = vessel(:,4);
N = length(r);
r_smooth = movmean(r,5);
diff_r_smooth = abs(diff(r_smooth));
ID = find(diff_r_smooth == min(diff_r_smooth));
canulla_r   = r_smooth(ID);
figure; subplot(2,1,1); hold on;
plot(r,'k','LineWidth',3);
plot(r_smooth,'--r','LineWidth',3);
plot(1:N,canulla_r.*ones(N,1),':c');
subplot(2,1,2); hold on;
plot(diff_r_smooth,'LineWidth',3);
alpha       = (scaling)./(2.*canulla_r);
end

