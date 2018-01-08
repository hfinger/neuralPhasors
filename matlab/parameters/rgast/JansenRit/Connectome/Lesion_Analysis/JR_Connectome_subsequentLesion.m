%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script that simulates behavior of connectome via a network of delay-coupled % 
% neural mass models with subsequent lesioning of critical nodes              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear paramsAll;
clear params;

%% set parameters of grid job class

params.Gridjob.runLocal = false; % if true, run all jobs on the pc the script is started from
params.Gridjob.requiremf = 3000; % required memory in mega byte
params.Gridjob.wc_host = ''; % names of PC's that will be excluded from cluster
params.Gridjob.jobname = 'JR_Connectome_subsequentLesion_'; % name of the grid job
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true; % use different random seeds if true
params.Gridjob.combParallel = false; % if true, grid search parameters need same number of entries
params.Gridjob.requiredThreads = '3'; % threads required per job

%% set parameters of specific class

params.JansenRitConnectome.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome.k = 3; %global connection strength scaling
params.JansenRitConnectome.v = 2; % velocity [m/s]
params.JansenRitConnectome.tMax = 605; %max simulation time [seconds]
params.JansenRitConnectome.dt = 0.0005; % simulation step size [seconds]
params.JansenRitConnectome.sampling = 2; % sampling every x steps
params.JansenRitConnectome.noiseVar = 22; % variance of noise used to drive neural masses
params.JansenRitConnectome.noiseMu = 220; % mean of noise used to drive neural masses
params.JansenRitConnectome.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectome.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome.subInp = [0,1,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome.d = 304; %initial interval to remove [seconds]
params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectome.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectome.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectome.FCMeasure = {{'Coherence'}}; % FC measure to use
params.JansenRitConnectome.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectome.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
params.JansenRitConnectome.gammaParams = true; % if true, use 30 Hz oscillators, else 10 Hz
params.JansenRitConnectome.nLesions = 5; % number of subsequent lesions (first run will be no lesions)
params.JansenRitConnectome.lesioningMode = 'minmax'; % how to choose the node to lesion (random, minmax or correlate)
params.JansenRitConnectome.targetF = 'driver'; % determines, whether coherence should be evaluated at network or driving frequency

%% look for pairs of nodes with specified path lengths between them and use them to position extrinsic drivers

% get Connectome and paths between nodes
p = params.JansenRitConnectome.p;
[C,D] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 4; % used to binarize C - number of inputs each node can have
maxPathLength = 2; % maximum number of nodes a path is allowed to have not counting the starting node
nodes = 1:33; % nodes for which to build paths
[~,idx] = sort(C,2); % sort C after input strength
C_binarized = zeros(size(C));
for i=1:length(nodes)
    C_binarized(i,idx(i,end-threshold+1:end)) = 1; % set n = threshold strongest connections to 1, rest 0
end
[~,~,~,~,paths,~] = findpaths(C_binarized, nodes, maxPathLength, 1); % build paths for each starting node

% set C and D
params.JansenRitConnectome.C = C_binarized;
params.JansenRitConnectome.D = D;

% find pairs of nodes which are connected by at least one indirect
% connection over a single third node and have no direct connection
node_pairs = cell(0);
idx = triu(ones(length(nodes)),1); % disregard reversed node combinations
k = 1;
for i=1:length(nodes)
    for j=1:length(nodes)
        if idx(i,j) == 1
            idx_col = paths(1,:) == i;
            if sum(paths(2,idx_col) == j) > 0
                continue
            end
            if sum(paths(3,idx_col) == j) > 0
                node_pairs{k} = [i,j];
                k = k+1;
            end
        end
    end
end

%% set driver parameters

% set driver params
params.JansenRitConnectome.drivPos = node_pairs; % indices of network nodes to be driven
params.JansenRitConnectome.drivScale = 0.1; % driver strength [V]
params.JansenRitConnectome.drivPO = [0:0.25:1]*2*pi; % phase offset of drivers
params.JansenRitConnectome.drivFreq = 28.5; % frequency of driver [Hz]
params.JansenRitConnectome.drivDur = 605; % duration of driver [s]

%% start grid job

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);