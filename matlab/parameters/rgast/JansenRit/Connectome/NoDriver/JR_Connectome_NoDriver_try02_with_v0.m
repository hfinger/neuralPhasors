clear paramsAll;

clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.Gridjob.jobname = 'JR_Connectome_NoDriver_try02_with_v0';
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:59:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectomePaper.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectomePaper.k = 30; %global connection strength scaling
params.JansenRitConnectomePaper.v = 3.2; % velocity [m/s]
params.JansenRitConnectomePaper.tMax = 605; %max simulation time [seconds]
params.JansenRitConnectomePaper.dt = 0.0005; % simulation step size [seconds]
params.JansenRitConnectomePaper.sampling = 2; % sampling every x steps
params.JansenRitConnectomePaper.noiseVar = 22; % variance of noise used to drive neural masses
params.JansenRitConnectomePaper.noiseMu = 220; % mean of noise used to drive neural masses
params.JansenRitConnectomePaper.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectomePaper.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.subInp = [0,1,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.d = 124; %initial interval to remove [seconds]
params.JansenRitConnectomePaper.verbose = false; % if we want to print time steps to console

params.JansenRitConnectomePaper.drivFreq = 10; % frequency of driver [Hz]
params.JansenRitConnectomePaper.drivPos = [1,2]; % indices of network nodes to be driven
params.JansenRitConnectomePaper.drivPO = [0]; % phase offset of drivers
params.JansenRitConnectomePaper.drivScale = 0; % driver strength [V]
params.JansenRitConnectomePaper.drivDur = 0; % duration of driver [s]

params.JansenRitConnectomePaper.FCMeasure = {{'Coherence'}}; % FC measure to use
params.JansenRitConnectomePaper.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectomePaper.corrSimFC = true; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectomePaper.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectomePaper.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectomePaper.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectomePaper.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
params.JansenRitConnectomePaper.plasticity = false; % if true, use 30 Hz oscillators, else 10 Hz

params.JansenRitConnectomePaper.u0 = 6e-3; % was 0 in master thesis!!!!!!

p = params.JansenRitConnectomePaper.p;
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
params.JansenRitConnectomePaper.C = bsxfun(@rdivide,C,sum(C,2));
params.JansenRitConnectomePaper.D = D;

% figure(1);
% imagesc(params.JansenRitConnectomePaper.C);
% figure(2);
% imagesc(params.JansenRitConnectomePaper.D);
% colormap(flipud(colormap()));

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);