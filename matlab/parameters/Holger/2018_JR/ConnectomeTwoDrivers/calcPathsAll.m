clear all;
p=1;
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
C = bsxfun(@rdivide,C,sum(C,2));
nodes = 1:33; % nodes for which to build paths
pairs = nchoosek(nodes,2);
maxPathLength = 9; % maximum number of nodes a path is allowed to have not counting the starting node
minPathNum = 1; % minimum number of paths between two nodes
[node_pairs,shortest_path_length,paths_all] = getPaths(C,pairs(:,1),pairs(:,2),maxPathLength,minPathNum);

data = dataPaths();
path_results = fullfile(data.resultsdir, 'Holger/2018_JR/ConnectomeTwoDrivers');
save(fullfile( path_results, 'paths.mat'), 'node_pairs','shortest_path_length','paths_all');