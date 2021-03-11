Initial inputs (outputs from DGTal create_distance_map and analyze_graph functions):
		  points_file.txt = txt file containing (x, y, z) coordinates of nodes, and edge connections with (x, y, z) coordinates of edge points 
                        (output from analyze_graph function, DGTal Github, https://github.com/DGtal-team/DGtal)
		  dmap.nrrd = distance map file (output from create_distance_map function, DGTal Github, https://github.com/DGtal-team/DGtal)

1) [arcs, nodes] = formatData('points_file.txt', 'dmap.nrrd');

% This program converts the information from the text and NRRD files into structures in Matlab:
% nodes: [# nodes x 5] matrix with columns [node ID | x | y | z | degree]
% arcs: {1 x # arcs} cell array
% 	-> arcs{1,i} = [(# edge points + 3) x 4] matrix
%		-> arcs{1,i}(1,:) 	  = 	  [  node1 ID |   node2 ID |       0      |        0      ]
%		-> arcs{1,i}(2,:) 	  = 	  [  X_node1  |   Y_node1  |    Z_node1   |   Rad_node1   ]
%		-> arcs{1,i}(3:(end-1),:) = [X_edgept_j | Y_edgept_j |  Z_edgept_j  |  Rad_edgept_j ] for 1 <= j <= # edge points
%		-> arcs{1,i}(end,:)   	=  	[  X_node2  |   Y_node2  |    Z_node2   |   Rad_node2   ]
% Employs function 'nrrdread.m'. Jeff Mather (2020). NRRD Format File Reader (https://www.mathworks.com/matlabcentral/fileexchange/34653-nrrd-format-file-reader), MATLAB Central File Exchange. Retrieved January 1, 2019.

2) [root_nodeID]=plotSlicerData(arcs, nodes)

% This program plots the nodes and edge points from the arcs file and labels the nodes with their node IDs.
% From this, we identify the following:
% root_nodeID = nodeID of the root node of the tree, the node that starts the MPA (taken to be node with min z coordinate)
% edges_to_remove = list of [node1 node2] pairs for false branches which need to be removed (visually identify)
% nodes_to_remove = list of nodeIDs which need to be removed from the tree (visually identify)
%			(Do not include nodes which are now degree 2 but still need to be points in the tree. We will eliminate these later.)

3) [path, arcs, nodes]=correctionEngine(arcs_old, nodes_old, root_node)

% This program works to correct 4 error types: false branches (FB), double edges (DE), duplicate points (DP), and small cycles (SC).
% FB: Requires inputs for edges and nodes which need to be removed from the graph.
% DE: If 2 edges exist from node A to node B, removes the larger of the two.
% DP: Removes any duplicates of edge points in the tree.
% SC: Removes edges of short cycles that do not show up in the shortest path from the root to a node.
% In addition, this program also merges all nodes which, after correction, now have degree=2.
% Outputs include updated nodes and arcs arrays, plus a new structure called "path".
% path: {1 x # nodes} cell array
% 	-> path{1,i} = [1 x # nodes from root] vector. 
%		-> path{1,i}(1) = root_nodeID
%		-> path{1,i}(end) = i^th nodeID ( nodes(i,1) )
%		-> path{1,i}(2:(end-1) = list of nodeIDs which comprise the shortest path from the root node to i^th node
% Employs function 'dijkstra.m'. Joseph Kirk (2020). Dijkstra's Shortest Path Algorithm (https://www.mathworks.com/matlabcentral/fileexchange/12850-dijkstra-s-shortest-path-algorithm), MATLAB Central File Exchange. Retrieved January 1, 2019.

4) [orientation, newNetwork, connectivity, nodeDetails] = directedGraphInfo(arcs, nodes, path)

% This program computes the orientation of each vessel in arcs, uses that to generate the newNetwork matrix of vessel info, and uses that to generate the connectivity matrix of node info. We also compute the nodeDetails of the tree, containing information about any trifurcations, quadfurcations, etc.
% orientation: [# arcs x 1] vector
%		->orientation(i) = 1 if node1=arcs{1,i}(1,1) comes before node2=arcs{1,i}(1,2) in the path
%				   = -1 if node2=arcs{1,i}(1,2) comes before node1=arcs{1,i}(1,1) in the path
% newNetwork: [# arcs x 6] matrix with columns [vessel ID | start nodeID | end nodeID | vessel radius | vessel length | radii std dev]
% 		->vessel ID of a vessel represented by arcs{1,i} is precisely i.
%		->vessel radius = interquartile mean of radius estimates at each node and edge point along that vessel
%		->vessel length = sum of Euclidean distances between consecutive points along the vessel
%		->radii std dev = standard deviation of the radius measurements at each node and edge point along that vessel
% connection: [# nodes x (4 + maxNumDaughters)] with columns [nodeID | out degree | in degree | daughter vessel IDs | parent vessel ID]
% 		->daughter vessel IDs = maxNumDaughters columns containing daughter vessel IDs (usually 2-4 columns)
% maxNumDaughters: maximum number of daughters present at any single junction in the network
% maxNumParents: maximum number of parents present at any single junction in the network (should always be 1 for networks which are trees)
% Employs the functions 'edge_orientation.m', 'networkGenerator.m', 'IQM.m', and 'connection.m'.


