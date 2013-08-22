clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'connEnv';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.ConnectomeEnvelope.inFileRates = 'connRate/4SimResult.mat';
params.ConnectomeEnvelope.t_rm = 20;
params.ConnectomeEnvelope.sigBandpass{1}.Fst1 = 1.5;
params.ConnectomeEnvelope.sigBandpass{1}.Fp1 = 2;
params.ConnectomeEnvelope.sigBandpass{1}.Fp2 = 6;
params.ConnectomeEnvelope.sigBandpass{1}.Fst2 = 6.5;
params.ConnectomeEnvelope.sigBandpass{1}.Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass{1}.Ap = 1;
params.ConnectomeEnvelope.sigBandpass{1}.Ast2 = 40;
params.ConnectomeEnvelope.sigBandpass{2}.Fst1 = 10;
params.ConnectomeEnvelope.sigBandpass{2}.Fp1 = 10.5;
params.ConnectomeEnvelope.sigBandpass{2}.Fp2 = 21.5;
params.ConnectomeEnvelope.sigBandpass{2}.Fst2 = 22;
params.ConnectomeEnvelope.sigBandpass{2}.Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass{2}.Ap = 1;
params.ConnectomeEnvelope.sigBandpass{2}.Ast2 = 40;
params.ConnectomeEnvelope.envLowpass.Fp = 0.5;
params.ConnectomeEnvelope.envLowpass.Fst = 1;
params.ConnectomeEnvelope.envLowpass.Ap = 1;
params.ConnectomeEnvelope.envLowpass.Ast = 40;
params.ConnectomeEnvelope.filtermethod = 'butter'; %or equiripple
params.ConnectomeEnvelope.outFilenames = 'fastFC';
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


