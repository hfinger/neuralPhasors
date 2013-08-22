clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'Kuramoto';
params.Gridjob.initRandStreamWithSeed = 1;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.GridToyKuramoto.network.weightedNetwork=false;
params.GridToyKuramoto.network.binaryNetwork=true;
params.GridToyKuramoto.network.N=5;
params.GridToyKuramoto.network.k_intraCluster = 0.9;
params.GridToyKuramoto.network.delay_intraClust = {2,4}; %in ms
params.GridToyKuramoto.network.numCluster = 3;
params.GridToyKuramoto.network.k_interCluster = 0.1;
params.GridToyKuramoto.network.delay_interClust = 5; %in ms

params.GridToyKuramoto.sim.k=12;
params.GridToyKuramoto.sim.f=40;
params.GridToyKuramoto.sim.v=10;
params.GridToyKuramoto.sim.t_max=10;
params.GridToyKuramoto.sim.dt=0.0001;
params.GridToyKuramoto.sim.sampling=10;
params.GridToyKuramoto.sim.sig_n=0;
params.GridToyKuramoto.sim.d=0;
params.GridToyKuramoto.sim.verbose=true;
params.GridToyKuramoto.sim.approx=false;
params.GridToyKuramoto.sim.invertSin=false;

params.GridToyKuramoto.env.t_rm=4;

params.GridToyKuramoto.env.sigBandpass(1).Fst1 = 5.5;
params.GridToyKuramoto.env.sigBandpass(1).Fp1 = 6;
params.GridToyKuramoto.env.sigBandpass(1).Fp2 = 22;
params.GridToyKuramoto.env.sigBandpass(1).Fst2 = 22.5;
params.GridToyKuramoto.env.sigBandpass(1).Ast1 = 40;
params.GridToyKuramoto.env.sigBandpass(1).Ap = 1;
params.GridToyKuramoto.env.sigBandpass(1).Ast2 = 40;

params.GridToyKuramoto.env.sigBandpass(2).Fst1 = 29.5;
params.GridToyKuramoto.env.sigBandpass(2).Fp1 = 30;
params.GridToyKuramoto.env.sigBandpass(2).Fp2 = 48;
params.GridToyKuramoto.env.sigBandpass(2).Fst2 = 48.5;
params.GridToyKuramoto.env.sigBandpass(2).Ast1 = 40;
params.GridToyKuramoto.env.sigBandpass(2).Ap = 1;
params.GridToyKuramoto.env.sigBandpass(2).Ast2 = 40;

params.GridToyKuramoto.env.envLowpass.Fp = 0.5;
params.GridToyKuramoto.env.envLowpass.Fst = 1;
params.GridToyKuramoto.env.envLowpass.Ap = 1;
params.GridToyKuramoto.env.envLowpass.Ast = 40;
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


