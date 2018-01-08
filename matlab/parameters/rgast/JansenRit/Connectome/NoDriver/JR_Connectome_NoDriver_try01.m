clear paramsAll;

clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.Gridjob.jobname = 'JR_Connectome_NoDriver_try01';
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:59:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectome2018.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome2018.k = 30; %global connection strength scaling
params.JansenRitConnectome2018.v = 3.2; % velocity [m/s]
params.JansenRitConnectome2018.tMax = 605; %max simulation time [seconds]
params.JansenRitConnectome2018.dt = 0.0005; % simulation step size [seconds]
params.JansenRitConnectome2018.sampling = 2; % sampling every x steps
params.JansenRitConnectome2018.noiseVar = 22; % variance of noise used to drive neural masses
params.JansenRitConnectome2018.noiseMu = 220; % mean of noise used to drive neural masses
params.JansenRitConnectome2018.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectome2018.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome2018.subInp = [0,1,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectome2018.d = 124; %initial interval to remove [seconds]
params.JansenRitConnectome2018.verbose = false; % if we want to print time steps to console

params.JansenRitConnectome2018.drivFreq = 10; % frequency of driver [Hz]
params.JansenRitConnectome2018.drivPos = [1,2]; % indices of network nodes to be driven
params.JansenRitConnectome2018.drivPO = [0]; % phase offset of drivers
params.JansenRitConnectome2018.drivScale = 0; % driver strength [V]
params.JansenRitConnectome2018.drivDur = 0; % duration of driver [s]

params.JansenRitConnectome2018.FCMeasure = {{'Coherence'}}; % FC measure to use
params.JansenRitConnectome2018.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectome2018.corrSimFC = true; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome2018.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectome2018.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectome2018.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectome2018.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
params.JansenRitConnectome2018.plasticity = false; % if true, use 30 Hz oscillators, else 10 Hz

p = params.JansenRitConnectome2018.p;
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
params.JansenRitConnectome2018.C = bsxfun(@rdivide,C,sum(C,2));
params.JansenRitConnectome2018.D = D;

% figure(1);
% imagesc(params.JansenRitConnectome2018.C);
% figure(2);
% imagesc(params.JansenRitConnectome2018.D);
% colormap(flipud(colormap()));

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);