clear paramsAll;

eeg.subj.Ids = [1:4 6:10];
eeg.subj.Avg = true;
eeg.day.Ids = [];
eeg.day.Avg = true;
eeg.cond.Ids = 5:6;
eeg.cond.Avg = true;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG_onlyMatchWithAicohcoh_Fisher2';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {10,11,12,13,14};
params.ConnectomeEnvelopeReduce.compareSC = true;
params.ConnectomeEnvelopeReduce.calcFisher = 2;

params.ConnectomeEnvelopeReduce.eeg = eeg;

params.ConnectomeEnvelopeReduce.sim.subjId.Ids = [];
params.ConnectomeEnvelopeReduce.sim.subjId.Avg = true;
% params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
% params.ConnectomeEnvelopeReduce.sim.k.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_onlyMatchWithAicohcoh_Fisher2';
params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder = [];%cellfun(@(x) ['db' num2str(x)],params.ConnectomeEnvelopeReduce.eegDatabase,'UniformOutput',false);

params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no';
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.calcCompareAicohWithCoh = true;

params.ConnectomeEnvelopeReduce.reloadConnFC = [];%'../CompareWithEEG_ConnFC/ConnFC.mat';
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [];
params.ConnectomeEnvelopeReduce.plotPermTests = false;
params.ConnectomeEnvelopeReduce.excludePlotFieldsRegexp = 'overFreq';
params.ConnectomeEnvelopeReduce.deletePlotFolder = false;
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll(2));
start(gridjobs);

