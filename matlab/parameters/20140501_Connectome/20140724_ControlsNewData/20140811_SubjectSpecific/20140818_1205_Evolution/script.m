clear paramsAll;
clear paramsAllRange;

subjects = [1:4];

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 2; % 0=datasimu from Arnaud, 1=Bastian, 2=BastianNew
params.ConnectomeSim.subjId = num2cell(subjects);

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = 0.02;
params.ConnectomeSim.normRow = 1;
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
params.ConnectomeSim.k=500;
params.ConnectomeSim.v=1.8;
params.ConnectomeSim.delay=1.5; %typically 0.3-0.5 ms up to 2 ms
params.ConnectomeSim.t_max={190,190};
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=false;
params.ConnectomeSim.statsRemoveInitialT = 10;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll(1).params = params;

paramsAllRange(1).params.ConnectomeSim.k = [200 1400 false true false];
paramsAllRange(1).params.ConnectomeSim.v = [0.3 14 false true false];
paramsAllRange(1).params.ConnectomeSim.delay = [0 10 false true false];
paramsAllRange(1).params.ConnectomeSim.homotopic = [0 0.5 false true false];

[variableParams, paramComb] = getVariableParams(paramsAll(1).params,false);

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeEnvelope';
params.Gridjob.initRandStreamWithSeed = 1;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.ConnectomeEnvelope.inFileRates = cellfun(@(x) ['ConnectomeSim/' num2str(x) 'SimResult.mat'], num2cell(1:size(paramComb,2)),'UniformOutput',false);
params.ConnectomeEnvelope.source_t_start = -Inf;
params.ConnectomeEnvelope.source_t_end = Inf;
params.ConnectomeEnvelope.saveSamples_t_start = [20 55];
params.ConnectomeEnvelope.saveSamples_t_end = [25 60];
params.ConnectomeEnvelope.env_t_start = 10;
params.ConnectomeEnvelope.env_t_end = Inf;
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

params.ConnectomeEnvelope.outFilenames = 'ConnectomeEnvelope';
params.ConnectomeEnvelope.saveSourceRate = false;
params.ConnectomeEnvelope.saveSourcePhase = false;
params.ConnectomeEnvelope.saveSourceBP = false;
params.ConnectomeEnvelope.saveEnvSig = false;
params.ConnectomeEnvelope.saveEnvLP = false;
paramsAll(2).params = params;
paramsAllRange(2).params = struct();
% paramsAllRange(2).params.ConnectomeEnvelope.env_t_start = [4 6 false true];

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'CompareWithEEG';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = 4;

params.ConnectomeEnvelopeReduce.eeg.subj.Ids = subjects;
params.ConnectomeEnvelopeReduce.eeg.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg.cond.Avg = true;

params.ConnectomeEnvelopeReduce.sim.t_max.Ids = [];
params.ConnectomeEnvelopeReduce.sim.t_max.Avg = true;

params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'match'; %reduce two subject dims to one dim
params.ConnectomeEnvelopeReduce.results.subj.Ids = [];
params.ConnectomeEnvelopeReduce.results.subj.Avg = true; % average over remaining subj dim

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithRS';
params.ConnectomeEnvelopeReduce.doPlot = false;
params.ConnectomeEnvelopeReduce.deleteEnvelopeWhenDone = true;
params.ConnectomeEnvelopeReduce.deleteSimWhenDone = true;
paramsAll(3).params = params;
paramsAllRange(3).params = struct();

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'Evolution';
params.Gridjob.continue = true;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = false;
params.Evolution.paramsAll = paramsAll;
params.Evolution.paramsAllRange = paramsAllRange;
params.Evolution.numGenerations = Inf;
params.Evolution.populationSize = 10;
params.Evolution.numNewIndividualsPerGen = 5;
params.Evolution.startGenId = 1;
params.Evolution.mutationInitial = 0.4;
params.Evolution.mutationHalfIter = 3;
params.Evolution.initFromOneIndividual = true;
params.Evolution.fitnessFcn = @(individualPath) max(0,subsref(load(fullfile(individualPath,'CompareWithRS','compareSimExp1.mat')), substruct('.','perFreqAvg','.','coh','.','rho')));
params.Evolution.doReplaceLeastFitParent = true;
evolutionParams{1} = params;

gridjobs = Gridjob(evolutionParams);
start(gridjobs);
