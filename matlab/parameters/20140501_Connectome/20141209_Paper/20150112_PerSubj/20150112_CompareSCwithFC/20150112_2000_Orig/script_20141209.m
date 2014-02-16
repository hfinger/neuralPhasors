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
params.ConnectomeEnvelopeReduce.eegDatabase = 20;

params.ConnectomeEnvelopeReduce.compareSC = true;
params.ConnectomeEnvelopeReduce.compareSC_db = 4;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_ConnFC';
params.ConnectomeEnvelopeReduce.onlyCollectConnFC = true;

paramsAll{1} = params;


eeg.subj.Ids = [1:4 6:10 11:13 15 17:20];
eeg.subj.Avg = false;
eeg.day.Ids = [];
eeg.day.Avg = true;
eeg.cond.Ids = 5:6;
eeg.cond.Avg = true;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {24,25};
params.ConnectomeEnvelopeReduce.compareSC = true;
params.ConnectomeEnvelopeReduce.compareSC_db = 4;

params.ConnectomeEnvelopeReduce.eeg = eeg;

params.ConnectomeEnvelopeReduce.sim.subjId.Ids = 1:17;
params.ConnectomeEnvelopeReduce.sim.subjId.Avg = false;
% params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
% params.ConnectomeEnvelopeReduce.sim.k.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG';
params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder = [];%cellfun(@(x) ['db' num2str(x)],params.ConnectomeEnvelopeReduce.eegDatabase,'UniformOutput',false);
params.ConnectomeEnvelopeReduce.savePreEvaluationData = true;

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

