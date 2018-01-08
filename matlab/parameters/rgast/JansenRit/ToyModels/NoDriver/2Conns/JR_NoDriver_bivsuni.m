clear paramsAll;

clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.Gridjob.jobname = 'JR_NoDriver_bivsuni';
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:59:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectome.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome.k = 10; %global connection strength scaling
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
params.JansenRitConnectome.plasticity = false;

% connectivity matrix
C = zeros(3);
C(1,3) = 1;
C(3,1) = 1;
C(1,2) = 1;
C(2,3) = 1;

% delay matrix
d_indirect = 100; % distance between outer nodes for indirect conn
d_direct = [0:200]; % distance between outer nodes for direct conn

for d=1:length(d_direct)
    D = zeros(3);
    D(1,2) = 0.5*d_indirect;
    D(2,3) = 0.5*d_indirect;
    D(1,3) = d_direct(d);
    D = D + D';
    Delays{d} = D;
end

params.JansenRitConnectome.C = C; % connectivity matrix
params.JansenRitConnectome.D = Delays; % distance matrix

params.JansenRitConnectome.drivPos = [1,3]; % indices of network nodes to be driven
params.JansenRitConnectome.drivScale = 0; % driver strength [V]
params.JansenRitConnectome.drivPO = 0; % phase offset of drivers
params.JansenRitConnectome.drivDur = 0; % duration of driver [s]

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);