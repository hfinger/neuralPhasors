%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script that simulates behavior of connectome via a network of delay-coupled % 
% neural mass models with a single lesion of each network node                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear paramsAll;
clear params;

%% set parameters of grid job class

params.Gridjob.runLocal = false; % if true, run all jobs on the pc the script is started from
params.Gridjob.wc_host = ''; % names of PC's that will be excluded from cluster
params.Gridjob.jobname = 'JR_Connectome_singleLesion_re'; % name of the grid job
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true; % use different random seeds if true
params.Gridjob.combParallel = true; % if true, grid search parameters need same number of entries
params.Gridjob.requiredThreads = '3'; % threads required per job
params.Gridjob.requiremf = 3000; % required memory in mega byte

%% set parameters of specific class

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
params.JansenRitConnectome.FCMeasure = {{'Coherence'},{'Coherence'},{'Coherence'},{'Coherence'}}; % FC measure to use
params.JansenRitConnectome.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectome.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
params.JansenRitConnectome.plasticity = false; % if true, use 30 Hz oscillators, else 10 Hz
params.JansenRitConnectome.lesionNode = num2cell([23,6,17,26]); % index of node to lesion in network

%% look for pairs of nodes with specified path lengths between them and use them to position extrinsic drivers

p = params.JansenRitConnectome.p;
[C,D] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1; % used to binarize C - number of inputs each node can have
nodes = 1:33; % nodes for which to build paths
C(C < threshold) = 0;
C = bsxfun(@rdivide,C,sum(C,2));

% set C and D
params.JansenRitConnectome.C = C;
params.JansenRitConnectome.D = D;

% find all paths and shortest path lengths between all pairs of nodes
% only consider non-directly connected pairs of nodes which connect via at
% most maxPathLength-1 other nodes and via at least 2 paths
pairs = nchoosek(nodes,2);
maxPathLength = 5; % maximum number of nodes a path is allowed to have not counting the starting node
minPathNum = 2; % minimum number of paths between two nodes
[node_pairs,shortest_path_length,paths_all] = getPaths(C,pairs(:,1),pairs(:,2),maxPathLength,minPathNum);

%% set driver parameters

% select random subset of node pairs for each path length
unique_path_lengths = unique(shortest_path_length);
n_node_pairs = 25;
driv_positions = cell(1,n_node_pairs*length(unique_path_lengths));
for i=1:length(unique_path_lengths)
    node_pairs_tmp = node_pairs(shortest_path_length == unique_path_lengths(i),:);
    idx = randperm(length(node_pairs_tmp),n_node_pairs);
    for j=1:n_node_pairs
        driv_positions{(i-1)*n_node_pairs+j} = node_pairs_tmp(idx(j),:);
    end
end

% set driver params
params.JansenRitConnectome.drivPos = {[4,13],[4,13],[4,13],[3,31]}; % indices of network nodes to be driven
params.JansenRitConnectome.drivScale = 20; % driver strength [mV]
params.JansenRitConnectome.drivPO = [0:0.0625:1]*2*pi; % phase offset of drivers
params.JansenRitConnectome.drivFreq = 11; % frequency of driver [Hz]
params.JansenRitConnectome.drivDur = 605; % duration of driver [s]

%% start grid job

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);