paths = dataPaths();
lf = load(fullfile(paths.databases,'SC_Bastian','surfaces','projectionMat','caAll_projectionMat.mat'));
lf.projMatRoi = mean(cell2mat(permute(lf.projMatRoi,[1 3 2])),3);
lf.projMatRoiPerVertex = mean(cell2mat(permute(lf.projMatRoiPerVertex,[1 3 2])),3);

clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 4; % 4 --> 17 Controls
params.ConnectomeSim.subjId = [1:4 6:10 11:13 15 17:20];

shufflePermutations = load('../../permuteTmp66.mat');
params.ConnectomeSim.shuffleSC = true;
params.ConnectomeSim.shuffleD = true;
params.ConnectomeSim.shufflePermutations = shufflePermutations.permuteTmp;

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = num2cell(0:0.02:0.22);
params.ConnectomeSim.normRow = 1;

params.ConnectomeSim.model = 'SAR'; % 'kuramoto' or 'rate'
params.ConnectomeSim.normStd = false;
params.ConnectomeSim.k = num2cell([0.4:0.05:0.95]);
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{1} = params;

[variableParams, paramComb] = getVariableParams(params,false);


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG_ConnFCLf';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = {'ConnectomeEnvelopeLfPerROI','ConnectomeEnvelopeLfPerVertex'};
params.ConnectomeEnvelopeReduce.eegDatabase = 26;

params.ConnectomeEnvelopeReduce.onlyFCsim = true;
params.ConnectomeEnvelopeReduce.FCsimApplyLf = {lf.projMatRoi, lf.projMatRoiPerVertex};
params.ConnectomeEnvelopeReduce.outDirectory = {'CompareWithEEG_ConnFCLfPerROI','CompareWithEEG_ConnFCLfPerVertex'};
params.ConnectomeEnvelopeReduce.onlyCollectConnFC = true;
params.ConnectomeEnvelopeReduce.useEnvFreqAsParamVar = false;

paramsAll{2} = params;

%%
eeg.subj.Ids = [1:4 6:10 11:13 15 17:20];
eeg.subj.Avg = true;
eeg.day.Ids = 1;
eeg.day.Avg = false;
eeg.cond.Ids = 5:6;
eeg.cond.Avg = true;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEGLf';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = 26;
params.ConnectomeEnvelopeReduce.compareSC = false;
params.ConnectomeEnvelopeReduce.onlyFCsim = false;
params.ConnectomeEnvelopeReduce.useEnvFreqAsParamVar = true;

params.ConnectomeEnvelopeReduce.eeg = eeg;

params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = [];
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEGLf';
params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder = [];

params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no';
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.calcCompareAicohWithCoh = true;

params.ConnectomeEnvelopeReduce.reloadConnFC = {'../CompareWithEEG_ConnFCLfPerROI/ConnFC.mat','../CompareWithEEG_ConnFCLfPerVertex/ConnFC.mat'};
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [];
params.ConnectomeEnvelopeReduce.plotPermTests = false;
params.ConnectomeEnvelopeReduce.excludePlotFieldsRegexp = 'overFreq';
params.ConnectomeEnvelopeReduce.deletePlotFolder = false;
params.ConnectomeEnvelopeReduce.doPlot = false;
paramsAll{3} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);