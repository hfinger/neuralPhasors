clear paramsAll;

% for the revision running stimulation on all nodes we change the
% simulation time parameters as follows:
% tMax = 90 (was 605)
% drivDur = 90 (was 605)
% initSampRem = 30 (was 124)


clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '!(*ramsauer*)';%'!(*saturn*|*c01grid*|*ramsauer*|*kuma*)';
params.Gridjob.jobname = 'Connectome';
params.Gridjob.runOnlyJobIds = [];
params.Gridjob.runOnlyUnfinishedJobIds = true;
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.combParallel = false;
params.Gridjob.walltime = '00:59:59';
params.Gridjob.requiredThreads = '3';

params.JansenRitConnectomePaper.p = 1; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectomePaper.k = 14; %num2cell([3, 10]); %30; %num2cell(round(22:2:34)); %global connection strength scaling
params.JansenRitConnectomePaper.v = 2.6; %num2cell(2.4:0.1:3.2); %3.2; % velocity [m/s]
params.JansenRitConnectomePaper.tMax = 70; %max simulation time [seconds]
params.JansenRitConnectomePaper.dt = 0.001; % simulation step size [seconds]
params.JansenRitConnectomePaper.sampling = 1; % sampling every x steps
params.JansenRitConnectomePaper.use_rk4 = true;
params.JansenRitConnectomePaper.noiseVar = 0; % variance of noise used to drive neural masses (22 was mentioned in msc thesis)
params.JansenRitConnectomePaper.noiseMu = [120, 320]; %num2cell(round(20:20:140)); %220; % mean of noise used to drive neural masses (220 was mentioned in msc thesis)
params.JansenRitConnectomePaper.runningAvg = false; % if true, substract running average from input to neural masses
params.JansenRitConnectomePaper.netInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.subInp = [1,0,0]; % scales the input to [pyramidal, excitatory, inhibitory] neurons from connectome
params.JansenRitConnectomePaper.initSampRem = 10; %initial interval to remove [seconds]
params.JansenRitConnectomePaper.verbose = true; % if we want to print time steps to console

c_tmp = 135;
params.JansenRitConnectomePaper.cs = [c_tmp, c_tmp*0.8, c_tmp*0.25, c_tmp*0.25, 0]; % connectivity strength (can sometimes be interpreted as average synaptic contacts)

params.JansenRitConnectomePaper.fTarget = 'drivFreq'; % [Hz]

params.JansenRitConnectomePaper.FCMeasure = {{'Coherence'}}; % FC measure to use
params.JansenRitConnectomePaper.fullFC = true; %whether to store full coherence matrix or only coherence of driven region
params.JansenRitConnectomePaper.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectomePaper.nWindows = 1; %number of timewindows over which to calculate mean coherence
params.JansenRitConnectomePaper.nBins = 32; % number of bins to use to discretize phase signal (only for MI and SE)
params.JansenRitConnectomePaper.storeY = true; %if true, store raw signal on simResults
params.JansenRitConnectomePaper.filterSig = true; % if true, raw PSPs will be bandpass-filtered before coherence is calculated

params.JansenRitConnectomePaper.u0 = 6e-3; % membrane voltage for which 50 % of maximum mean firing rate is observed [V].. was 0 in master thesis
params.JansenRitConnectomePaper.He = 3.25e-3; % Average synaptic gain for excitatory synapses [V]
params.JansenRitConnectomePaper.Hi = 22e-3; % Average synaptic gain for inhibitory synapses [V]
params.JansenRitConnectomePaper.Te = 10e-3; % Average time constant for excitatory signal transfer (synaptic delays,..) [s]
params.JansenRitConnectomePaper.Ti = 20e-3; % Average time constant for inhibitory signal transfer (synaptic delays,..) [s]

% Te and Ti correspond to 1/a and 1/b in the Jansen-Rit Paper.

params.JansenRitConnectomePaper.e0 = 2.5; % determines maximum mean firing rate [1/s]
params.JansenRitConnectomePaper.r = 560; % steepness of the sigmoid transfer function from mean psp to mean firing rate [1/V]

params.JansenRitConnectomePaper.washout = 2000;

params.JansenRitConnectomePaper.calcFC_nwin1 = false;
params.JansenRitConnectomePaper.subtract_S0 = false;
params.JansenRitConnectomePaper.use_moran = false;
params.JansenRitConnectomePaper.use_out_psp = false;
params.JansenRitConnectomePaper.use_sigm_as_out = false;
params.JansenRitConnectomePaper.use_inpP_as_out = true;
params.JansenRitConnectomePaper.use_sigm_y0_as_out = false;
params.JansenRitConnectomePaper.calcCohWithDriver = true;
params.JansenRitConnectomePaper.calcKuramotoOrderParam = true;
      
p = params.JansenRitConnectomePaper.p;
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
params.JansenRitConnectomePaper.C = bsxfun(@rdivide,C,sum(C,2));
params.JansenRitConnectomePaper.D = D;


% set driver params
%idx = randperm(length(C),2);
idx = 1;%length(C);

params.JansenRitConnectomePaper.drivPosVarMatrix = idx'; % indices of network nodes to be driven
params.JansenRitConnectomePaper.drivPos = []; % indices of network nodes to be driven
params.JansenRitConnectomePaper.drivScale = 0.25; % driver strength [mV]
params.JansenRitConnectomePaper.drivPO = 0; % phase offset of drivers
params.JansenRitConnectomePaper.drivFreq = 31; % frequency of driver [Hz]
params.JansenRitConnectomePaper.drivDur = 90; % duration of driver [s]

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);