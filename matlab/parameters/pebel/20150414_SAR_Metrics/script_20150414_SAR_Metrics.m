clear paramsAll; 

% global representation of parameter values
glob.a = num2cell(0.0:0.1:0.7);
glob.b = num2cell(1:1:2);
glob.c = num2cell(1:1:4);
glob.d = num2cell(0.4:0.05:0.95);

%%

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeMetrics';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeMetrics.sparse = glob.a;
params.ConnectomeMetrics.heuristics = glob.b;
params.ConnectomeMetrics.hScale = glob.c;
params.ConnectomeMetrics.graphH0 = 0; 
paramsAll{1} = params;

%%

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 5;                                          
params.ConnectomeSim.loadSp = glob.a;                         
params.ConnectomeSim.loadHeur = glob.b;                         
params.ConnectomeSim.loadHscale = glob.c;

params.ConnectomeSim.subjId = [1:4 6:10 11:13 15 17:20];

% these params are redundant, corresp. functions are implemented in ConnectomeMetrics:
params.ConnectomeSim.normRowBeforeHomotopic = 0;
params.ConnectomeSim.homotopic = 0;                                        
params.ConnectomeSim.normRow = 0;                                          

params.ConnectomeSim.model = 'SAR';                                         % 'kuramoto' or 'rate'
params.ConnectomeSim.normStd = true;                                        % true: correlation output, false: covariance output
params.ConnectomeSim.k = glob.d;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{2} = params;

%%

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

% parameters to command calculation of local prediction errors
params.ConnectomeEnvelopeReduce.calcSquaredDist = false; % calcSquaredDistAvg 
params.ConnectomeEnvelopeReduce.sim = struct();
% see: function [tensor, outSpec] = filterTensor(tensor, filterSpec, inSpec )

paramsAll{3} = params;

%%

eeg.subj.Ids = [1:4 6:10 11:13 15 17:20];                                   % same as params.ConnectomeSim.subjId
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
params.ConnectomeEnvelopeReduce.compareSC = false;                          % compare FC directly with SC
params.ConnectomeEnvelopeReduce.compareSC_db = 5;                           % value 5 to compare with homotopic dataset
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

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'ConnectomeAnalysis';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;
params.Gridjob.runOnlyJobIds = [5];

params.ConnectomeAnalysis.loadSp = glob.a;                         
params.ConnectomeAnalysis.loadHeur = glob.b;                         
params.ConnectomeAnalysis.loadHscale = glob.c;
params.ConnectomeAnalysis.loadKscale = glob.d;
paramsAll{5} = params;

% write a script for post-simulation analysis that plots metrics as fct of parameters
% re-assess global and local prediction errors
% get plots, organize results in reasonable folder structure

%%
clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);