clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithSeed = 2;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 1; % 0=datasimu from Arnaud, 1=SC from Bastian
params.ConnectomeSim.subjId = 2;%{2,3,4,6,7,8,9};

params.ConnectomeSim.normRow = 1;
params.ConnectomeSim.model = 'rate'; % 'kuramoto' or 'rate'
params.ConnectomeSim.useNetworkFokkerPlanck = false;

%params specific for Kuramoto:
params.ConnectomeSim.approx=false;
params.ConnectomeSim.invertSin=false;
params.ConnectomeSim.f=40;
params.ConnectomeSim.startState = [];

%params specific for rate model:
params.ConnectomeSim.tau={10,20,40};

%params for all models:
params.ConnectomeSim.k={0.1, 1, 10, 100};
params.ConnectomeSim.v={1, 3, 10};
params.ConnectomeSim.t_max=80;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n={0, 0.2};
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 20;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{1} = params;

% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
% params.Gridjob.jobname = 'boldSig';
% params.Gridjob.initRandStreamWithSeed = 12345;
% params.Gridjob.continue = false;
% params.ConnectomeFCeval.inFileRates = {'connRate/1.mat','connRate/2.mat','connRate/3.mat'};
% params.ConnectomeFCeval.t_rm = 20;
% params.ConnectomeFCeval.outFilenames = 'boldSig';
% paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


