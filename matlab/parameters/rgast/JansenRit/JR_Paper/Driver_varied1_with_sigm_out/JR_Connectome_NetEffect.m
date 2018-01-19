clear paramsAll;

clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.Gridjob.jobname = 'JR_Connectome_NetEffect';
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:59:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectome.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome.k = 30; %global connection strength scaling
params.JansenRitConnectome.v = 3.2; % velocity [m/s]
params.JansenRitConnectome.tMax = 605; %max simulation time [seconds]
params.JansenRitConnectome.dt = 0.0005; % simulation step size [seconds]
params.JansenRitConnectome.sampling = 2; % sampling every x steps
params.JansenRitConnectome.noiseVar = 22; % variance of noise used to drive neural masses
params.JansenRitConnectome.noiseMu = 220; % mean of noise used to drive neural masses
params.JansenRitConnectome.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectome.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome.subInp = [0,1,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome.d = 124; %initial interval to remove [seconds]
params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectome.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectome.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectome.FCMeasure = {{'Coherence'}}; % FC measure to use
params.JansenRitConnectome.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectome.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
params.JansenRitConnectome.plasticity = false; % if true, add synaptic plasticity mechanism to pyramidal cell synapses

%% look for pairs of nodes with specified path lengths between them and use them to position extrinsic drivers

p = params.JansenRitConnectome.p;
[C,D] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1; % used to binarize C - number of inputs each node can have
maxPathLength = 5; % maximum number of nodes a path is allowed to have not counting the starting node
nodes = 1:33; % nodes for which to build paths
C(C < threshold) = 0;
C = bsxfun(@rdivide,C,sum(C,2));

% set C and D
params.JansenRitConnectome.C = C;
params.JansenRitConnectome.D = D;

% find paths between each pair of nodes (non-symmetric). store paths, 
% shortest path lengths and pairs
paths_all = cell(0);
node_pairs = cell(0);
shortest_path_length = cell(0);
n = 1;
for i=1:length(nodes)
    [~,~,~,~,paths,~] = findpaths(C > 0, i, maxPathLength, 1); % find all paths
    for j=1:length(nodes)
        idx_col = find(sum(paths == j,1) > 0); % find all paths from node i to node j
        if j == i || C(i,j) > 0 || C(j,i) > 0 || isempty(idx_col) % disregard direct and self connections
            continue
        else
            idx_row = zeros(1,length(idx_col));
            for k=1:length(idx_col) 
                idx_row(k) = find(paths(:,idx_col(k)) == j); % find position of node j on path
                paths(idx_row(k)+1:end,idx_col(k)) = 0; % set all following nodes to 0
            end
            paths_all{1,n} = paths(:,idx_col);
            shortest_path_length{n} = min(idx_row) - 2;
            node_pairs{n,1} = [i,j];
            n = n+1;
        end
    end
end

node_pairs = cell2mat(node_pairs);

% check for bi-directional double entries
for i=1:size(node_pairs,1)
    paths = unique(paths_all{i}','rows')';
    if size(paths,2) < 2
        node_pairs(i,:) = 0; % disregard all node pairs with less than two connections
    end
    idx = all(node_pairs == repmat(fliplr(node_pairs(i,:)),size(node_pairs,1),1),2);
    node_pairs(idx,:) = 0; % disregard flipped node pairs
end
idx = node_pairs(:,1) > 0;
node_pairs = node_pairs(idx,:);
shortest_path_length = {shortest_path_length{idx}};
shortest_path_length = cell2mat(shortest_path_length);

%% set driver parameters

% select random subset of node pairs for each path length
unique_path_lengths = unique(shortest_path_length);
n_node_pairs = 1;
driv_positions = cell(1,n_node_pairs*length(unique_path_lengths));
for i=1:length(unique_path_lengths)
    node_pairs_tmp = node_pairs(shortest_path_length == unique_path_lengths(i),:);
    idx = randperm(length(node_pairs_tmp),n_node_pairs);
    for j=1:n_node_pairs
        driv_positions{(i-1)*n_node_pairs+j} = node_pairs_tmp(idx(j),:);
    end
end

% set driver params
params.JansenRitConnectome.drivPos = driv_positions; % indices of network nodes to be driven
params.JansenRitConnectome.drivScale = num2cell([10,20,30,40,50]); % driver strength [V]
params.JansenRitConnectome.drivPO = num2cell([0:0.0625:1]*2*pi); % phase offset of drivers
params.JansenRitConnectome.drivFreq = {11,11,11,11}; % frequency of driver [Hz]
params.JansenRitConnectome.drivDur = 605; % duration of driver [s]

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);