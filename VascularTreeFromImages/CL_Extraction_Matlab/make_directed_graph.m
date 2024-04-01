function make_directed_graph(conn,terminal)
%% Takes simple connectivity and displays a directed graph
% ASSUMES SIMPLE BIFURCATIONS ONLY
% MJC 12/21/17
conn = conn;
M = size(conn,1);
T = length(terminal);
starts = zeros(1+2*M,1);
ends   = zeros(1+2*M,1);
weights = zeros(1+2*M,1);
strahler = zeros(max(max(conn))-1,1);
starts(1) = 1;
ends(1) = 2;
weights(1) = 2;
for i=1:M
    starts(2*i)     = conn(i,1);
    ends(2*i)       = conn(i,2);
    weights(2*i)    = conn(i,2);
    starts(2*i+1)   = conn(i,1);
    ends(2*i+1)     = conn(i,3);
    weights(2*i+1)  = conn(i,3);
    if any(terminal == conn(i,2))
        strahler(2*i) = 1;
    end
    if any(terminal == conn(i,3))
        strahler(2*i+1) = 1;
    end
end
weights = weights-2;
%Implement Strahler Order Number?
for i=2*M+1:-2:2 %%Walk backwards to inlet 
    d1 = strahler(i);
    d2 = strahler(i-1);
    if d1 == d2
        strahler(i-2) = d1+1;
    else
        strahler(i-2) = max(d1,d2)
    end
end
G = digraph(starts,ends);%,weights);
% figure(1234);clf; 
% title('Directed Graph using vessel number')
% h = plot(G,'EdgeLabel',G.Edges.Weight);
% h.NodeColor = 'red';
% h.LineWidth = 3;
% h.MarkerSize = 8;
% h.EdgeAlpha = 0.8;
% set(gca,'FontSize',30);

% G = digraph(starts,ends,strahler);
% figure(12345);clf;
% title('Directed Graph with Strahler Number');
% h = plot(G,'EdgeLabel',G.Edges.Weight);
% h.NodeColor = 'red';
% h.LineWidth = 3;
% h.MarkerSize = 8;
% h.EdgeAlpha = 0.8;
% set(gca,'FontSize',30);

end




% % Custom colors from parula colormap and custom sizes between 12 and 20 pts
% labelColors = parula(length(names));
% labelSizes  = randi([12 20], length(names), 1);
% p.NodeLabel = {};
% % Custom labels
% hold on
% for i=1:length(names)
%     text(p.XData(i), p.YData(i), names(i), 'Color', labelColors(i, :), 'FontSize', labelSizes(i));
% end
% hold off