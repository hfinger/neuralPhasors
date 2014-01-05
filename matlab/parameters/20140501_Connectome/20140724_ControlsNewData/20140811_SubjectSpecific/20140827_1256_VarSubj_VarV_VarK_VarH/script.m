clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 2; % 0=datasimu from Arnaud, 1=Bastian, 2=BastianNew
params.ConnectomeSim.subjId = num2cell([1:4 6:10]);

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = num2cell(0:0.05:0.1);
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
params.ConnectomeSim.k=num2cell(400:100:1000);
params.ConnectomeSim.v=num2cell(1:0.5:4);
params.ConnectomeSim.delay=1.5; %typically 0.3-0.5 ms up to 2 ms
params.ConnectomeSim.t_max={100,100};
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 10;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{1} = params;

[variableParams, paramComb] = getVariableParams(paramsAll{1},false);

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
params.ConnectomeEnvelope.env_t_start = 20;
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
paramsAll{2} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'CompareWithEEG';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {2,3,4};

params.ConnectomeEnvelopeReduce.eeg{1}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{1}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{1}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{1}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{1}.prepost.Ids = 1;
params.ConnectomeEnvelopeReduce.eeg{1}.prepost.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{1}.task.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{1}.task.Avg = true;

params.ConnectomeEnvelopeReduce.eeg{2}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{2}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{2}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{2}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{2}.prepost.Ids = 1;
params.ConnectomeEnvelopeReduce.eeg{2}.prepost.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{2}.task.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{2}.task.Avg = true;

params.ConnectomeEnvelopeReduce.eeg{3}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{3}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{3}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{3}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{3}.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg{3}.cond.Avg = true;

params.ConnectomeEnvelopeReduce.sim.t_max.Ids = [];
params.ConnectomeEnvelopeReduce.sim.t_max.Avg = true;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_allSim';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'match'; % 'no' or 'match' or 'nonmatch'
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [3 4 2 1];
params.ConnectomeEnvelopeReduce.permutePlotDimsPermTest = [2 3 1];
params.ConnectomeEnvelopeReduce.plotPermTests = true;

paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'CompareWithEEG';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {2,3,4};

params.ConnectomeEnvelopeReduce.eeg{1}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{1}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{1}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{1}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{1}.prepost.Ids = 1;
params.ConnectomeEnvelopeReduce.eeg{1}.prepost.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{1}.task.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{1}.task.Avg = true;

params.ConnectomeEnvelopeReduce.eeg{2}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{2}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{2}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{2}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{2}.prepost.Ids = 1;
params.ConnectomeEnvelopeReduce.eeg{2}.prepost.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{2}.task.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{2}.task.Avg = true;

params.ConnectomeEnvelopeReduce.eeg{3}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{3}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{3}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{3}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{3}.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg{3}.cond.Avg = true;

%% TODO: change ids
params.ConnectomeEnvelopeReduce.sim.t_max.Ids = [];
params.ConnectomeEnvelopeReduce.sim.t_max.Avg = true;
params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = 1;
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;
params.ConnectomeEnvelopeReduce.sim.v.Ids = [];
params.ConnectomeEnvelopeReduce.sim.v.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_allSubjPairs';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no'; % 'no' or 'match' or 'nonmatch'
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [];
params.ConnectomeEnvelopeReduce.plotPermTests = true;

paramsAll{4} = params;


clear params;
params.Gridjob.runLocal = true;
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

params.ConnectomeEnvelopeReduce.eeg.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg.cond.Avg = true;

params.ConnectomeEnvelopeReduce.sim.t_max.Ids = [];
params.ConnectomeEnvelopeReduce.sim.t_max.Avg = true;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_new';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'match'; % 'no' or 'match' or 'nonmatch'
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [3 4 2 1];
params.ConnectomeEnvelopeReduce.permutePlotDimsPermTest = [2 3 1];
params.ConnectomeEnvelopeReduce.plotPermTests = false;

paramsAll{5} = params;

clear params;
params.Gridjob.runLocal = true;
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

params.ConnectomeEnvelopeReduce.eeg.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg.cond.Avg = true;

params.ConnectomeEnvelopeReduce.sim.t_max.Ids = [];
params.ConnectomeEnvelopeReduce.sim.t_max.Avg = true;
params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = 1;
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = 5;
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;
params.ConnectomeEnvelopeReduce.sim.v.Ids = 5;
params.ConnectomeEnvelopeReduce.sim.v.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_newBest';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no'; % 'no' or 'match' or 'nonmatch'
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [];
params.ConnectomeEnvelopeReduce.permutePlotDimsPermTest = [];
params.ConnectomeEnvelopeReduce.plotPermTests = false;

paramsAll{5} = params;


clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer';
params.Gridjob.jobname = 'CompareWithEEG';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {2,3,4,6,6,2,3,4,6,6};

params.ConnectomeEnvelopeReduce.eeg{1}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{1}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{1}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{1}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{1}.prepost.Ids = 1;
params.ConnectomeEnvelopeReduce.eeg{1}.prepost.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{1}.task.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{1}.task.Avg = true;

params.ConnectomeEnvelopeReduce.eeg{2}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{2}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{2}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{2}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{2}.prepost.Ids = 1;
params.ConnectomeEnvelopeReduce.eeg{2}.prepost.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{2}.task.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{2}.task.Avg = true;

params.ConnectomeEnvelopeReduce.eeg{3}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{3}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{3}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{3}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{3}.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg{3}.cond.Avg = true;

% real RS lcmv:
params.ConnectomeEnvelopeReduce.eeg{4}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{4}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{4}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{4}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{4}.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg{4}.cond.Avg = true;

% task:
params.ConnectomeEnvelopeReduce.eeg{5}.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg{5}.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg{5}.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg{5}.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg{5}.cond.Ids = 1:2;
params.ConnectomeEnvelopeReduce.eeg{5}.cond.Avg = true;

params.ConnectomeEnvelopeReduce.eeg(6:10) = params.ConnectomeEnvelopeReduce.eeg(1:5);

params.ConnectomeEnvelopeReduce.sim.t_max.Ids = [];
params.ConnectomeEnvelopeReduce.sim.t_max.Avg = true;
params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = 1;
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;
params.ConnectomeEnvelopeReduce.sim.v.Ids = [];
params.ConnectomeEnvelopeReduce.sim.v.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_newAll';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = [repmat({'match'},1,5) repmat({'no'},1,5)];
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [repmat({[2 3 1]},1,5) repmat({[1 2 3 4]},1,5)];
params.ConnectomeEnvelopeReduce.plotPermTests = false;

paramsAll{6} = params;

clear params;
gridjobs = Gridjob(paramsAll(6));
start(gridjobs);

