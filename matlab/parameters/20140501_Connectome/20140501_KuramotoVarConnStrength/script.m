clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'Kuramoto';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.Gridjob.remoteStart = true;

params.ConnectomeSim.model = 'kuramoto';
params.ConnectomeSim.k={0.01 , 0.1, 0.3, 0.9, 5};
params.ConnectomeSim.f=60;
params.ConnectomeSim.v=10;
params.ConnectomeSim.t_max=500;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=1.25;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.approx=false;
params.ConnectomeSim.outFilenames = 'results';

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


