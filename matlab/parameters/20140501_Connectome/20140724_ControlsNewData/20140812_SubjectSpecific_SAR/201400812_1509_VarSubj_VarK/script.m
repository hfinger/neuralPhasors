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
params.ConnectomeSim.homotopic = 0.02;
params.ConnectomeSim.normRow = 1;

params.ConnectomeSim.model = 'SAR'; % 'kuramoto' or 'rate'
params.ConnectomeSim.k = num2cell(0.05:0.05:0.95);
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{1} = params;

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

params.ConnectomeEnvelopeReduce.outDirectory = {'CompareWithPreTask','CompareWithRS'};
params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'match';
params.ConnectomeEnvelopeReduce.results.subj.Ids = [];
params.ConnectomeEnvelopeReduce.results.subj.Avg = true;

paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll{2});
start(gridjobs);

