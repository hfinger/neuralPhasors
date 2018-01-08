clear paramsAll;

clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'JR_2Driver_k';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';

params.JansenRitConnectome.p = num2cell([2]); %defines what kind of p norm to use for normalization of structural connectivity
params.JansenRitConnectome.k = num2cell([0.03,0.1,0.3,1,3,10,30]); %global connection strength scaling
params.JansenRitConnectome.v = num2cell([4]); % velocity [m/s]
params.JansenRitConnectome.tMax = 125; %max simulation time [seconds]
params.JansenRitConnectome.dt = 0.0001; % simulation step size [seconds]
params.JansenRitConnectome.sampling = 10; % sampling every x steps
params.JansenRitConnectome.snr = num2cell([0]); % amount of noise
params.JansenRitConnectome.d = 5; %initial interval to remove [seconds]
params.JansenRitConnectome.verbose = false; % if we want to print time steps to console
params.JansenRitConnectome.drivFreq = num2cell([29.5]); %driving frequencies
params.JansenRitConnectome.drivPos = [4,31]; % network node to drive
params.JansenRitConnectome.drivRange = [20,15]; % variance of the gaussian centered around DrivPos determining the strength of stimulation
params.JansenRitConnectome.drivPO = [0:0.25:1.75]*pi; % phase offset of drivers
params.JansenRitConnectome.drivScale = num2cell([0.02]); % amplitude/strength of drivers
params.JansenRitConnectome.drivStart = 1; % timepoint at which to start driving [seconds]
params.JansenRitConnectome.drivDur = 125; % driving duration [seconds]
params.JansenRitConnectome.corrSimFC = false; % if true, compare simulated coherence to empirical coherence
params.JansenRitConnectome.fullCoherence = false; %whether to store full coherence matrix or only coherence of driven region

paramsAll{1} = params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);