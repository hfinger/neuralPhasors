clear paramsAll;

clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'JR_2Driver_v';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;
params.Gridjob.runOnlyJobIds = [];

params.JansenRitConnectome.p = 2; %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome.k = 8; %global connection strength scaling
params.JansenRitConnectome.v = num2cell([1,2,4,8]); % velocity [m/s]
params.JansenRitConnectome.tMax = 15; %max simulation time [seconds]
params.JansenRitConnectome.dt = 0.0001; % simulation step size [seconds]
params.JansenRitConnectome.sampling = 10; % sampling every x steps
params.JansenRitConnectome.snr = 0; % amount of noise
params.JansenRitConnectome.d = 5; %initial interval to remove [seconds]
params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
params.JansenRitConnectome.drivFreq = 29.5; %driving frequencies
params.JansenRitConnectome.drivPos = [4,31]; % network node to drive
params.JansenRitConnectome.drivRange = [20,15]; % variance of the gaussian centered around DrivPos determining the strength of stimulation
params.JansenRitConnectome.drivPO = [0:0.5:1.]*pi; % phase offset of drivers
params.JansenRitConnectome.drivScale = 0.02]); % amplitude/strength of drivers
params.JansenRitConnectome.drivStart = 1; % timepoint at which to start driving [seconds]
params.JansenRitConnectome.drivDur = 15; % driving duration [seconds]
params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome.fullCoherence = false; %whether to store full coherence matrix or only coherence of driven region

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);