clear paramsAll;

clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.Gridjob.jobname = 'JR_DelayEffects_test3_';
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:59:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectome.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome.k = 20; %global connection strength scaling
params.JansenRitConnectome.v = 2; % velocity [m/s]
params.JansenRitConnectome.tMax = 105; %max simulation time [seconds]
params.JansenRitConnectome.dt = 0.0005; % simulation step size [seconds]
params.JansenRitConnectome.sampling = 2; % sampling every x steps
params.JansenRitConnectome.noiseVar = 22; % variance of noise used to drive neural masses
params.JansenRitConnectome.noiseMu = 220; % mean of noise used to drive neural masses
params.JansenRitConnectome.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectome.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome.subInp = [0,1,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome.d = 5; %initial interval to remove [seconds]
params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectome.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectome.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectome.FCMeasure = {{'Coherence','SE'}}; % FC measure to use
params.JansenRitConnectome.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectome.filterSig = false; % if true, raw PSPs will be bandpass-filtered before coherence is calculated

% connectivity matrix
C = zeros(3);
Connectivities{1} = C;

C(1,3) = 1;
C(3,1) = 1;
Connectivities{2} = C;

C(1,2) = 1;
C(2,1) = 1;
Connectivities{3} = C;

C(2,3) = 1;
C(3,2) = 1;
Connectivities{4} = C;

C(1,3) = 0;
C(3,1) = 0;
Connectivities{5} = C;

C(2,3) = 0;
C(3,2) = 0;
Connectivities{6} = C;

C(1,2) = 0;
C(2,1) = 0;
C(2,3) = 1;
C(3,2) = 1;
Connectivities{7} = C;

C(1,3) = 1;
C(3,1) = 1;
Connectivities{8} = C;

% delay matrix
d12 = [10:10:160]; % distance between nodes at the ends of chain
D = zeros(3);
D(1,3) = 100;
D(2,3) = 40;
for d=1:length(d12)
    D_tmp = D;
    D_tmp(1,2) = d12(d);
    D_tmp = D_tmp + D_tmp';
    Delays{d} = D_tmp;
end

params.JansenRitConnectome.C = Connectivities{2}; % connectivity matrix
params.JansenRitConnectome.D = Delays; % distance matrix
params.JansenRitConnectome.drivPos = [1,3]; % indices of network nodes to be driven
params.JansenRitConnectome.drivScale = 0.1; % driver strength [V]
params.JansenRitConnectome.drivPO = num2cell([0:0.125:1.9375]*pi); % phase offset of drivers {0.5*pi, 1.5*pi};%
params.JansenRitConnectome.drivFreq = 9; % frequency of driver [Hz]
params.JansenRitConnectome.drivStart = 0; % starting point of driver [s]
params.JansenRitConnectome.drivDur = 605; % duration of driver [s]

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);