% Computes some connectivity information.. Make this a part of
% the above function in the future.
function [vessel_info] = update_details_cl(vessel_info, iteration, numvessels, conn)
% len = 0;
% disp(iteration)
vessel = vessel_info{iteration, 2};
if isempty(vessel)
    vessel_info{iteration,1} = 'EMPTY VESSEL';
    return
end
%% Added on 2/19/17
% Check to see if the vessel is going from inlet to outlet or visa versa.
% diff_ves_z = diff(vessel(:,3));
% if mean(diff_ves_z) < 0
%     vessel = flipud(vessel);
% else
% end
%%
x = vessel(:,1);
y = vessel(:,2);
z = vessel(:,3);
r = vessel(:,4);
N = length(r);
L = length(x);
if L < 10
    cl_leng = 0;
    for i=1:L-1
        cl_leng = cl_leng + sqrt( (x(i) - x(i+1)).^2 +  ...
            (y(i) - y(i+1)).^2 + (z(i) - z(i+1)).^2);
    end
    radii = mean(r);
    median_radii = median(r);
    r_std = std(r);
    inlet = r(1);
    outlet = r(end);
    taper = outlet./inlet;
else
    cl_leng = 0;
    for i=1:L-1
        cl_leng = cl_leng + sqrt( (x(i) - x(i+1)).^2 +  ...
            (y(i) - y(i+1)).^2 + (z(i) - z(i+1)).^2);
    end
    quarter_radius = round(size(r,1)./4);
    radii = mean(r(quarter_radius:3*quarter_radius));
    median_radii = median(r(quarter_radius:3*quarter_radius));
    r_std = std(r);
    upperB = round(0.9.*N); %Upper bound is last 10%
    lowerB = round(0.1.*N); %Lower bound is first 10%
    if lowerB > 0
        inlet = mean(r(1:lowerB));
        outlet = mean(r(upperB:end));
    else
        inlet = r(1);
        outlet = r(end);
    end
    taper = outlet./inlet;
end
figure(232);
plot(r,'o');
vessel_info{iteration,2} = vessel;
vessel_info{iteration,3} = cl_leng;
vessel_info{iteration,4} =  radii;
vessel_info{iteration,5} = median_radii;
vessel_info{iteration,6} = taper;
vessel_info{iteration,7} = inlet;
vessel_info{iteration,8} = outlet;
vessel_info{iteration,9} = r_std;
% Now find daughter and parent
flag1 = 0; % For Daughter
flag2 = 0; % For Parent
i = 2;
i_end = numvessels;

while flag1 + flag2 < 2
    if i == i_end
        if flag2 == 0
            
            vessel_info{iteration,10} = 'NO PARENT';
            %         elseif flag2 ==0
            %             vessel_details{1,7} = 'TERMINAL'
            
            break
        end
    end
    
    if flag1 == 0
        if strcmp(vessel_info{iteration,1}, conn{i,1})
            vessel_info{iteration,11} = strcat(conn{i,2}, ',', conn{i,3}, ...
                ',', conn{i,4}, ',', conn{i,5}, ',', conn{i,6}, ',', conn{i,7});
            flag1 = 1;
        end
    end
    
    
    if flag2 == 0
        if strcmp(vessel_info{iteration,1},conn{i,2}) || strcmp(vessel_info{iteration,1},conn{i,3}) ...
                || strcmp(vessel_info{iteration,1},conn{i,4}) || strcmp(vessel_info{iteration,1}, conn{i,5}) ...
                || strcmp(vessel_info{iteration,1},conn{i,6}) || strcmp(vessel_info{iteration,1},conn{i,7})
            vessel_info{iteration,10} = strcat(conn{i,1});
            flag2 = 1;
        end
    end
    
    i = i+1;
end
end

