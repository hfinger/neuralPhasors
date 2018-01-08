clear paramsAll;

clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.Gridjob.jobname = 'JR_Connectome_NoDriver_final2018';
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:59:00';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectome.p = num2cell(ones(1,17)); %defines what kind of p norm to use for normalization of structural connectivity
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

params.JansenRitConnectome.drivFreq = 10; % frequency of driver [Hz]
params.JansenRitConnectome.drivPos = [1,2]; % indices of network nodes to be driven
params.JansenRitConnectome.drivPO = [0]; % phase offset of drivers
params.JansenRitConnectome.drivScale = 0; % driver strength [V]
params.JansenRitConnectome.drivDur = 0; % duration of driver [s]

params.JansenRitConnectome.FCMeasure = {{'Coherence'}}; % FC measure to use
params.JansenRitConnectome.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectome.corrSimFC = true; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectome.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectome.storeY = false; %if true, store raw signal on simResults
params.JansenRitConnectome.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated
params.JansenRitConnectome.plasticity = false; % if true, use 30 Hz oscillators, else 10 Hz

p = params.JansenRitConnectome.p{1};
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
params.JansenRitConnectome.C = bsxfun(@rdivide,C,sum(C,2));
params.JansenRitConnectome.D = D;

% figure(1);
% imagesc(params.JansenRitConnectome.C);
% figure(2);
% imagesc(params.JansenRitConnectome.D);
% colormap(flipud(colormap()));

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);