clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 8000;
% params.Gridjob.wc_host = '!(*ramsauer*)';
params.Gridjob.jobname = 'ConnectomeSimTest';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.requiredThreads = '2';
params.Gridjob.runOnlyJobIds = [25 26];
params.Gridjob.continue = true;
params.Gridjob.walltime = '6:00:00';

params.ConnectomeSim.dataset = 4; % 4 --> 19 Controls
params.ConnectomeSim.subjId = -1;

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = 0.12;
params.ConnectomeSim.normRow = 1;

shufflePermutations = load('../../permuteTmp66.mat');
params.ConnectomeSim.shuffleSC = false;
params.ConnectomeSim.shuffleD = false;
params.ConnectomeSim.shufflePermutations = [];%shufflePermutations.permuteTmp;

params.ConnectomeSim.model = 'kuramoto'; % 'kuramoto' or 'rate'
params.ConnectomeSim.useNetworkFokkerPlanck = false;

%params specific for Kuramoto:
params.ConnectomeSim.approx=false;
params.ConnectomeSim.invertSin=false;
params.ConnectomeSim.f=60;
params.ConnectomeSim.startState = [];

%params specific for rate model:
params.ConnectomeSim.tau=20;

%params for all models:
params.ConnectomeSim.k=700;
params.ConnectomeSim.v=1.7;
params.ConnectomeSim.delay=1.25; %typically 0.3-0.5 ms up to 2 ms
params.ConnectomeSim.t_max=490;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 10;
params.ConnectomeSim.outFilenames = 'ConnectomeSimTest';
params.ConnectomeSim.forceOverwrite = false;

ConnectomeEnvelopeParams.inFileRates = [];
ConnectomeEnvelopeParams.source_t_start = -Inf;
ConnectomeEnvelopeParams.source_t_end = Inf;
ConnectomeEnvelopeParams.saveSamples_t_start = [20 55];
ConnectomeEnvelopeParams.saveSamples_t_end = [25 60];
ConnectomeEnvelopeParams.env_t_start = 10;
ConnectomeEnvelopeParams.env_t_end = Inf;
ConnectomeEnvelopeParams.filtermethod = 'butter'; %or equiripple

ConnectomeEnvelopeParams.applyLeadField = true;
ConnectomeEnvelopeParams.applyLeadFieldMethod = 'perROI'; 

% paramm = 7; %Parmeter m see Tallon-Baudry & Bertrand, TICS, 1999
% centerFreq = 2.^(1:0.5:6);
% sigma_f = centerFreq / paramm;

% centerFreq = 1:30;
% sigma_f = 2*ones(size(centerFreq));
% sigma_f(1) = 0.1;
% sigma_f(2) = 1;

% centerFreq = [5 11 22];
% sigma_f = 2*ones(size(centerFreq));

centerFreq = [1 2 3 4 5 6 7 8 9 10 11 12];
sigma_f = 2*ones(size(centerFreq));
sigma_f(1) = 0.25;
sigma_f(2) = 1;

paramsAll{1} = params;

[variableParams, paramComb] = getVariableParams(paramsAll{1},false);






clear params;
gridjobs = Gridjob(paramsAll(1));
start(gridjobs);

