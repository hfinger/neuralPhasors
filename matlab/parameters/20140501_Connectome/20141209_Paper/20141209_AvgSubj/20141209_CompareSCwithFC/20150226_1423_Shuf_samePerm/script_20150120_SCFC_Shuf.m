clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG_ConnFC';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = 10;

params.ConnectomeEnvelopeReduce.compareSC = true;
params.ConnectomeEnvelopeReduce.compareSC_db = 4;

shufflePermutations = load('../../permuteTmp66.mat');
params.ConnectomeEnvelopeReduce.compareSC_shuffle = true;
params.ConnectomeEnvelopeReduce.compareSC_shufflePermutations = shufflePermutations.permuteTmp;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_ConnFC';
params.ConnectomeEnvelopeReduce.onlyCollectConnFC = true;

paramsAll{1} = params;


eeg.subj.Ids = [1:4 6:10 11:13 15 17:20];
eeg.subj.Avg = true;
eeg.day.Ids = 1;
eeg.day.Avg = false;
eeg.cond.Ids = 5:6;
eeg.cond.Avg = true;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG_onlyMatchWithAicohcoh';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {30,25,27,28,29};
params.ConnectomeEnvelopeReduce.compareSC = true;
params.ConnectomeEnvelopeReduce.compareSC_db = 4;

params.ConnectomeEnvelopeReduce.eeg = eeg;

params.ConnectomeEnvelopeReduce.sim.subjId.Ids = 1:17;
params.ConnectomeEnvelopeReduce.sim.subjId.Avg = true;
% params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
% params.ConnectomeEnvelopeReduce.sim.k.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG';
params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder = [];%cellfun(@(x) ['db' num2str(x)],params.ConnectomeEnvelopeReduce.eegDatabase,'UniformOutput',false);

params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no';
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.calcCompareAicohWithCoh = true;

params.ConnectomeEnvelopeReduce.reloadConnFC = '../CompareWithEEG_ConnFC/ConnFC.mat';
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [];
params.ConnectomeEnvelopeReduce.plotPermTests = false;
params.ConnectomeEnvelopeReduce.excludePlotFieldsRegexp = 'overFreq';
params.ConnectomeEnvelopeReduce.deletePlotFolder = false;
params.ConnectomeEnvelopeReduce.doPlot = false;
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll(2));
start(gridjobs);
