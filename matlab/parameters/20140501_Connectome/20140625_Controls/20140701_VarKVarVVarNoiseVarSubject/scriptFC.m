clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'connEnv';
params.Gridjob.initRandStreamWithSeed = 1;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.ConnectomeEnvelope.inFileRates = cellfun(@(x) ['ConnectomeSim/' num2str(x) 'SimResult.mat'], num2cell(1:126),'UniformOutput',false);
params.ConnectomeEnvelope.source_t_start = -Inf;
params.ConnectomeEnvelope.source_t_end = Inf;
params.ConnectomeEnvelope.saveSamples_t_start = [20 55];
params.ConnectomeEnvelope.saveSamples_t_end = [25 60];
params.ConnectomeEnvelope.env_t_start = 20;
params.ConnectomeEnvelope.env_t_end = 80;
params.ConnectomeEnvelope.filtermethod = 'butter'; %or equiripple

params.ConnectomeEnvelope.sigBandpass(1).Fst1 = 2.5;
params.ConnectomeEnvelope.sigBandpass(1).Fp1 = 3;
params.ConnectomeEnvelope.sigBandpass(1).Fp2 = 7;
params.ConnectomeEnvelope.sigBandpass(1).Fst2 = 7.5;
params.ConnectomeEnvelope.sigBandpass(1).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(1).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(1).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(2).Fst1 = 8.5;
params.ConnectomeEnvelope.sigBandpass(2).Fp1 = 9;
params.ConnectomeEnvelope.sigBandpass(2).Fp2 = 13;
params.ConnectomeEnvelope.sigBandpass(2).Fst2 = 13.5;
params.ConnectomeEnvelope.sigBandpass(2).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(2).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(2).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(3).Fst1 = 22.5;
params.ConnectomeEnvelope.sigBandpass(3).Fp1 = 23;
params.ConnectomeEnvelope.sigBandpass(3).Fp2 = 27;
params.ConnectomeEnvelope.sigBandpass(3).Fst2 = 27.5;
params.ConnectomeEnvelope.sigBandpass(3).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(3).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(3).Ast2 = 40;

params.ConnectomeEnvelope.envLowpass.Fp = 0.5;
params.ConnectomeEnvelope.envLowpass.Fst = 1;
params.ConnectomeEnvelope.envLowpass.Ap = 1;
params.ConnectomeEnvelope.envLowpass.Ast = 40;

params.ConnectomeEnvelope.outFilenames = 'connEnv';
params.ConnectomeEnvelope.saveSourceRate = true;
params.ConnectomeEnvelope.saveSourcePhase = true;
params.ConnectomeEnvelope.saveSourceBP = true;
params.ConnectomeEnvelope.saveEnvSig = true;
params.ConnectomeEnvelope.saveEnvLP = true;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


