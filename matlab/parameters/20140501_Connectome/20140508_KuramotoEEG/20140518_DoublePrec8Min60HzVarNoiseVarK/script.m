clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithSeed = 2;
params.Gridjob.continue = false;
params.Gridjob.exclusive = false;
params.Gridjob.requiredThreads = '3';
params.ConnectomeSim.normRow = true;
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
params.ConnectomeSim.k={100, 300, 1000, 3000};
params.ConnectomeSim.v=10;
params.ConnectomeSim.t_max=500;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n={0, 0.2, 0.8, 1.25};
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 20;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'boldSig';
params.Gridjob.initRandStreamWithSeed = 1;
params.Gridjob.exclusive = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.ConnectomeFCeval.inFileRates = cellfun(@(x) ['ConnectomeSim/' num2str(x) 'SimResult.mat'], num2cell(1:16),'UniformOutput',false);
params.ConnectomeFCeval.t_rm = 20;
params.ConnectomeFCeval.outFilenames = 'boldSig';
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll{2});
start(gridjobs);


