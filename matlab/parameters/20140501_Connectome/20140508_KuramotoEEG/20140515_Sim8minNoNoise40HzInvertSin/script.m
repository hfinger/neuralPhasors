clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'connRateInvert';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.ConnectomeSim.normRow = true;
params.ConnectomeSim.model = 'kuramoto';
params.ConnectomeSim.k={30, 100, 300};
params.ConnectomeSim.f=40;
params.ConnectomeSim.invertSin=true;
params.ConnectomeSim.v={2, 4, 8};
params.ConnectomeSim.t_max=100;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.approx=false;
params.ConnectomeSim.outFilenames = 'connRate';
params.ConnectomeSim.statsRemoveInitialT = 20;
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


