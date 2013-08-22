clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'connEnv';
params.Gridjob.initRandStreamWithSeed = 1;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '4';
params.ConnectomeEnvelope.inFileRates = cellfun(@(x) ['ConnectomeSim/' num2str(x) 'SimResult.mat'], num2cell(1:5),'UniformOutput',false);
params.ConnectomeEnvelope.source_t_start = -Inf;
params.ConnectomeEnvelope.source_t_end = Inf;
params.ConnectomeEnvelope.saveSamples_t_start = 1:120:500;
params.ConnectomeEnvelope.saveSamples_t_end = 3:120:500;
params.ConnectomeEnvelope.env_t_start = [20, 60];
params.ConnectomeEnvelope.env_t_end = [500, 500];
params.ConnectomeEnvelope.filtermethod = 'butter'; %or equiripple

params.ConnectomeEnvelope.sigBandpass(1).Fst1 = 1.5;
params.ConnectomeEnvelope.sigBandpass(1).Fp1 = 2;
params.ConnectomeEnvelope.sigBandpass(1).Fp2 = 6;
params.ConnectomeEnvelope.sigBandpass(1).Fst2 = 6.5;
params.ConnectomeEnvelope.sigBandpass(1).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(1).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(1).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(2).Fst1 = 3.5;
params.ConnectomeEnvelope.sigBandpass(2).Fp1 = 4;
params.ConnectomeEnvelope.sigBandpass(2).Fp2 = 8;
params.ConnectomeEnvelope.sigBandpass(2).Fst2 = 8.5;
params.ConnectomeEnvelope.sigBandpass(2).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(2).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(2).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(3).Fst1 = 5.5;
params.ConnectomeEnvelope.sigBandpass(3).Fp1 = 6;
params.ConnectomeEnvelope.sigBandpass(3).Fp2 = 10.5;
params.ConnectomeEnvelope.sigBandpass(3).Fst2 = 11;
params.ConnectomeEnvelope.sigBandpass(3).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(3).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(3).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(4).Fst1 = 7.5;
params.ConnectomeEnvelope.sigBandpass(4).Fp1 = 8;
params.ConnectomeEnvelope.sigBandpass(4).Fp2 = 13;
params.ConnectomeEnvelope.sigBandpass(4).Fst2 = 13.5;
params.ConnectomeEnvelope.sigBandpass(4).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(4).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(4).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(5).Fst1 = 10;
params.ConnectomeEnvelope.sigBandpass(5).Fp1 = 10.5;
params.ConnectomeEnvelope.sigBandpass(5).Fp2 = 21.5;
params.ConnectomeEnvelope.sigBandpass(5).Fst2 = 22;
params.ConnectomeEnvelope.sigBandpass(5).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(5).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(5).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(6).Fst1 = 12.5;
params.ConnectomeEnvelope.sigBandpass(6).Fp1 = 13;
params.ConnectomeEnvelope.sigBandpass(6).Fp2 = 30;
params.ConnectomeEnvelope.sigBandpass(6).Fst2 = 30.5;
params.ConnectomeEnvelope.sigBandpass(6).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(6).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(6).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(7).Fst1 = 21;
params.ConnectomeEnvelope.sigBandpass(7).Fp1 = 21.5;
params.ConnectomeEnvelope.sigBandpass(7).Fp2 = 39;
params.ConnectomeEnvelope.sigBandpass(7).Fst2 = 39.5;
params.ConnectomeEnvelope.sigBandpass(7).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(7).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(7).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(8).Fst1 = 29.5;
params.ConnectomeEnvelope.sigBandpass(8).Fp1 = 30;
params.ConnectomeEnvelope.sigBandpass(8).Fp2 = 48;
params.ConnectomeEnvelope.sigBandpass(8).Fst2 = 48.5;
params.ConnectomeEnvelope.sigBandpass(8).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(8).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(8).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(9).Fst1 = 38.5;
params.ConnectomeEnvelope.sigBandpass(9).Fp1 = 39;
params.ConnectomeEnvelope.sigBandpass(9).Fp2 = 66;
params.ConnectomeEnvelope.sigBandpass(9).Fst2 = 66.5;
params.ConnectomeEnvelope.sigBandpass(9).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(9).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(9).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(10).Fst1 = 51.5;
params.ConnectomeEnvelope.sigBandpass(10).Fp1 = 52;
params.ConnectomeEnvelope.sigBandpass(10).Fp2 = 80;
params.ConnectomeEnvelope.sigBandpass(10).Fst2 = 80.5;
params.ConnectomeEnvelope.sigBandpass(10).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(10).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(10).Ast2 = 40;

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


