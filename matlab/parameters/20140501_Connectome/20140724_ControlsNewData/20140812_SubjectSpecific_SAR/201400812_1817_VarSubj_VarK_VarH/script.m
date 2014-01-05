clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 2; % 0=datasimu from Arnaud, 1=Bastian, 2=BastianNew
params.ConnectomeSim.subjId = num2cell([1:4 6:10]);

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = num2cell(0:0.01:0.12);
params.ConnectomeSim.normRow = 1;

params.ConnectomeSim.model = 'SAR'; % 'kuramoto' or 'rate'
params.ConnectomeSim.k = num2cell([0.05:0.05:0.95 0.98]);
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{1} = params;

[variableParams, paramComb] = getVariableParams(params,false);

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
params.ConnectomeEnvelopeReduce.onlyFCsim = true;
params.ConnectomeEnvelopeReduce.eegDatabase = {3,4};

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
params.ConnectomeEnvelopeReduce.eeg{2}.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg{2}.cond.Avg = true;

params.ConnectomeEnvelopeReduce.outDirectory = {'matchingSubjParis/CompareWithPreTask','matchingSubjParis/CompareWithRS'};
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'match';
params.ConnectomeEnvelopeReduce.results.subj.Ids = [];
params.ConnectomeEnvelopeReduce.results.subj.Avg = true;

paramsAll{2} = params;


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
params.ConnectomeEnvelopeReduce.onlyFCsim = true;
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

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_allSim';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = [repmat({'match'},1,5) repmat({'no'},1,5)];
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = true;
params.ConnectomeEnvelopeReduce.permutePlotDims = [repmat({[2 3 1]},1,5) repmat({[1 2 3 4]},1,5)];
params.ConnectomeEnvelopeReduce.plotPermTests = false;

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
params.ConnectomeEnvelopeReduce.onlyFCsim = true;
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

params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = 11; % --> homotopic=0.1
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = 14;  % --> k=0.7
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_bestSim';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no'; % 'no' or 'match' or 'nonmatch'
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [];
params.ConnectomeEnvelopeReduce.plotPermTests = false;

paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'CompareWithEEG_linRegress';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.onlyFCsim = true;
params.ConnectomeEnvelopeReduce.eegDatabase = 5;

params.ConnectomeEnvelopeReduce.eeg.subj.Ids = [1:4 6:10];
params.ConnectomeEnvelopeReduce.eeg.subj.Avg = false;
params.ConnectomeEnvelopeReduce.eeg.day.Ids = [];
params.ConnectomeEnvelopeReduce.eeg.day.Avg = true;
params.ConnectomeEnvelopeReduce.eeg.cond.Ids = 5:6;
params.ConnectomeEnvelopeReduce.eeg.cond.Avg = true;

params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = 11; % --> homotopic=0.1
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = 14;  % --> k=0.7
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_linRegress';
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'match'; % 'no' or 'match' or 'nonmatch'
params.ConnectomeEnvelopeReduce.results = struct();
params.ConnectomeEnvelopeReduce.calcSquaredDist = true;

params.ConnectomeEnvelopeReduce.reloadConnFC = true;
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [];
params.ConnectomeEnvelopeReduce.plotPermTests = false;

paramsAll{5} = params;

clear params;
gridjobs = Gridjob(paramsAll(3));
start(gridjobs);