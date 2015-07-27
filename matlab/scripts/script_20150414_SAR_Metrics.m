% this script belongs into the 'parameters' directory, /src_pebel/matlab/parameters/20140501_Connectome/20141209_Paper/20141209_AvgSubj/20141209_SAR/20150120_1654_Orig_Params12x12

clear paramsAll; 

%%

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeMetrics';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeMetrics.sparse = num2cell(0.6:0.1:0.7);
params.ConnectomeMetrics.heuristics = num2cell(1:1:2);
params.ConnectomeMetrics.hScale = num2cell(1:1:2);
params.ConnectomeMetrics.graphH0 = false; 
paramsAll{1} = params;

%%

% [variableParams, paramComb] = getVariableParams(params,false);

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 5;                                           % load SC.mat and SCh.mat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.ConnectomeSim.subjId = [1:4 6:10 11:13 15 17:20];

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = 0;                                         % num2cell(0:0.02:0.22);
params.ConnectomeSim.normRow = 1;

params.ConnectomeSim.model = 'SAR'; % 'kuramoto' or 'rate'
params.ConnectomeSim.normStd = true;                                        % true: correlation output, false: covariance output
params.ConnectomeSim.k = num2cell(0.4:0.05:0.95);
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{2} = params;

%%

[variableParams, paramComb] = getVariableParams(params,false);

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

params.ConnectomeEnvelopeReduce.onlyFCsim = true;
params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_ConnFC';
params.ConnectomeEnvelopeReduce.onlyCollectConnFC = true;

paramsAll{3} = params;

%%

eeg.subj.Ids = [1:4 6:10 11:13 15 17:20]; %%%%%%%%%%%%%%%%%%%%%% same as params.ConnectomeSim.subjId
eeg.subj.Avg = true;
eeg.day.Ids = 1;  %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
eeg.day.Avg = false; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
params.Gridjob.runOnlyJobIds = [5];

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {24,25,27,28,29};
params.ConnectomeEnvelopeReduce.compareSC = false;
params.ConnectomeEnvelopeReduce.compareSC_db = 5; % value 5 to compare with homotopic dataset
params.ConnectomeEnvelopeReduce.onlyFCsim = true;

params.ConnectomeEnvelopeReduce.eeg = eeg;

params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = [];
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;

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
paramsAll{4} = params;

%%

% write a script for post-simulation analysis that plots metrics as fct of parameters
% re-assess global and local prediction errors
% get plots, organize results in reasonable folder structure

%%
clear params;
gridjobs = Gridjob(paramsAll(4));
start(gridjobs);