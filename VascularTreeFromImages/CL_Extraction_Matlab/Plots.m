load 'CLData.mat'


l_v_r = cell(length(vessel_details)-1,2);
for i = 2:155
    l = [linspace(0,vessel_details{i,3},height(vessel_details{i,2}))./10]';
    r = vessel_details{i,2}(1:end,4)./10;
    l_r = [l flip(r)];
    l_v_r{i-1,1} = vessel_details{i,1};
    l_v_r{i-1,2} = l_r;
end
% for i = 2:5
%     %l = linspace(0,vessel_details{i,3},height(vessel_details{i,2}));
%     l = [linspace(0,vessel_details{i,3},height(vessel_details{i,2}))./10]';
%     r = vessel_details{i,2}(1:end,4)./10;
%     % l_v_r = [l,flip(r)];
%     figure(i)
%     %subplot(2,1,1)
%     plot (l,flip(r),'.','MarkerSize',30)
%     xlabel('length');
%     ylabel('radius');
%     set(gca,'fontsize',20)
%     title(vessel_details{i,1})
% end
% 
% % for i = 2:5
% %     l = [linspace(0,vessel_details{i,3},height(vessel_details{i,2}))./10]';
% %     r = vessel_details{i,2}(1:end,4)./10;
% %     l_r = [l r]
% %     % l_v_r{i-1,1} = vessel_details{i,1};
% %     % l_v_r{i,2} = l_r
% % end
% 
% 
% for i = 2:5
%     %l = linspace(0,vessel_details{i,3},height(vessel_details{i,2}));
%     x = vessel_details{i,2}(1:end,3)./10;
%     y = vessel_details{i,2}(1:end,2)./10;
%     z = vessel_details{i,2}(1:end,1)./10;
%     figure(i)
%     subplot(2,1,2)
%     plot3(flip(x),flip(y),flip(z))
%     xlabel('x');
%     ylabel('y');
%     zlabel('z');
%     set(gca,'fontsize',20)
%     title(vessel_details{i,1})
% end

    % figure(vessel_details{2,1}, 2)
    % l = linspace(0,vessel_details{2,3},height(vessel_details{2,2}));
    % plot (l,vessel_details{2,2}(1:end,4))
    %height(vessel_details)