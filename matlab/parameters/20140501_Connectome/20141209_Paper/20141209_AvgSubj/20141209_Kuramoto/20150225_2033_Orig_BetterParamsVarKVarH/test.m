clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 8000;
params.Gridjob.wc_host = '!(*ramsauer*)';
params.Gridjob.jobname = 'testsdfsd';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.requiredThreads = '2';
params.Gridjob.runOnlyJobIds = [];
params.Gridjob.continue = true;
params.Gridjob.walltime = '6:00:00';

params.ConnectomeSim.dataset = 4; % 4 --> 19 Controls
params.ConnectomeSim.subjId = -1;

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = num2cell(0:0.02:0.22);
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
params.ConnectomeSim.k=num2cell(100:100:1200);
params.ConnectomeSim.v=1.7;%3.4;
params.ConnectomeSim.delay=1.25; %typically 0.3-0.5 ms up to 2 ms
params.ConnectomeSim.t_max=490;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 10;
params.ConnectomeSim.outFilenames = 'testsdfsd';
params.ConnectomeSim.forceOverwrite = false;

paramsAll{1} = params;

[variableParams, paramComb] = getVariableParams(paramsAll{1},false);



clear params;
gridjobs = Gridjob(paramsAll(1));
start(gridjobs);

