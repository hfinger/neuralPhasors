clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'Connectome';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.Gridjob.remoteStart = true;

params.ConnectomeSim.k=0.8;
params.ConnectomeSim.f=60;
params.ConnectomeSim.v=10;
params.ConnectomeSim.t_max=500;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n={0.1, 1.25, 3};
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.approx={false, true};
params.ConnectomeSim.t_rm = 20;
params.ConnectomeSim.outFilenames = 'results';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeFC';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.Gridjob.remoteStart = true;
params.ConnectomeFCeval.inFileRates = {'results/1.mat','results/2.mat','results/3.mat','results/4.mat','results/5.mat','results/6.mat'};
params.ConnectomeFCeval.t_rm = 20;
params.ConnectomeFCeval.outFilenames = 'boldSignal';
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll{2});
start(gridjobs);


