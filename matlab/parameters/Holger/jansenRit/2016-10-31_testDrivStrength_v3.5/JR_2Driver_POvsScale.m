clear paramsAll;

clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.wc_host = '';
params.Gridjob.jobname = 'JR_2Driver_testDrivStrength';
params.Gridjob.continue = false;
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = false;

params.JansenRitConnectome.p = 2; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome.k = 1; %global connection strength scaling
params.JansenRitConnectome.v = 3.5; % velocity [m/s]
params.JansenRitConnectome.tMax = 245; %max simulation time [seconds]
params.JansenRitConnectome.dt = 0.0001; % simulation step size [seconds]
params.JansenRitConnectome.sampling = 10; % sampling every x steps
params.JansenRitConnectome.snr = 0; % amount of noise
params.JansenRitConnectome.d = 5; %initial interval to remove [seconds]
params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
params.JansenRitConnectome.drivFreq = 29.5; %driving frequencies
params.JansenRitConnectome.drivPos = [4,31]; % network node to drive
params.JansenRitConnectome.drivRange = [20,15]; % variance of the gaussian centered around DrivPos determining the strength of stimulation
params.JansenRitConnectome.drivPO = [0:0.0625:1.9375]*pi; % phase offset of drivers
params.JansenRitConnectome.drivScale = num2cell(0:0.001:0.025); % amplitude/strength of drivers
params.JansenRitConnectome.drivStart = 1; % timepoint at which to start driving [seconds]
params.JansenRitConnectome.drivDur = 245; % driving duration [seconds]
params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome.fullCoherence = false; %whether to store full coherence matrix or only coherence of driven region
this.params.JansenRitConnectome.nWindows = 8; %number of timewindows over which to calculate mean coherence

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);