clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'connRate';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.ConnectomeSim.model = 'rate';
params.ConnectomeSim.k=0.8;
params.ConnectomeSim.tau=20;
params.ConnectomeSim.v=10;
params.ConnectomeSim.t_max=500;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n={0.1, 0.25, 0.5};
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.outFilenames = 'connRate';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'boldSig';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.ConnectomeFCeval.inFileRates = {'connRate/1.mat','connRate/2.mat','connRate/3.mat'};
params.ConnectomeFCeval.t_rm = 20;
params.ConnectomeFCeval.outFilenames = 'boldSig';
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


