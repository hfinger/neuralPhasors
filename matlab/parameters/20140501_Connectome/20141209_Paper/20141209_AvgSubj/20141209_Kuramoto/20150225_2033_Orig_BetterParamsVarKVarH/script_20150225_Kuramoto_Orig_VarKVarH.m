clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 8000;
params.Gridjob.wc_host = '!(*ramsauer*)';
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.requiredThreads = '2';
params.Gridjob.runOnlyJobIds = [];
params.Gridjob.continue = true;
params.Gridjob.walltime = '6:00:00';

params.ConnectomeSim.dataset = 4; % 4 --> 19 Controls
params.ConnectomeSim.subjId = -1;

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = num2cell(0:0.02:0.22);
params.ConnectomeSim.normRow = 1;

shufflePermutations = load('../../permuteTmp66.mat');
params.ConnectomeSim.shuffleSC = false;
params.ConnectomeSim.shuffleD = false;
params.ConnectomeSim.shufflePermutations = [];%shufflePermutations.permuteTmp;

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
params.ConnectomeSim.k=num2cell(100:100:1200);
params.ConnectomeSim.v=1.7;%3.4;
params.ConnectomeSim.delay=1.25; %typically 0.3-0.5 ms up to 2 ms
params.ConnectomeSim.t_max=490;
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 10;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
params.ConnectomeSim.forceOverwrite = false;

ConnectomeEnvelopeParams.inFileRates = [];
ConnectomeEnvelopeParams.source_t_start = -Inf;
ConnectomeEnvelopeParams.source_t_end = Inf;
ConnectomeEnvelopeParams.saveSamples_t_start = [20 55];
ConnectomeEnvelopeParams.saveSamples_t_end = [25 60];
ConnectomeEnvelopeParams.env_t_start = 10;
ConnectomeEnvelopeParams.env_t_end = Inf;
ConnectomeEnvelopeParams.filtermethod = 'butter'; %or equiripple

ConnectomeEnvelopeParams.applyLeadField = true;
ConnectomeEnvelopeParams.applyLeadFieldMethod = 'perROI'; 

% paramm = 7; %Parmeter m see Tallon-Baudry & Bertrand, TICS, 1999
% centerFreq = 2.^(1:0.5:6);
% sigma_f = centerFreq / paramm;

% centerFreq = 1:30;
% sigma_f = 2*ones(size(centerFreq));
% sigma_f(1) = 0.1;
% sigma_f(2) = 1;

% centerFreq = [5 11 22];
% sigma_f = 2*ones(size(centerFreq));

centerFreq = [1 2 3 4 5 6 7 8 9 10 11 12];
sigma_f = 2*ones(size(centerFreq));
sigma_f(1) = 0.25;
sigma_f(2) = 1;

for f=1:length(centerFreq)
  ConnectomeEnvelopeParams.sigBandpass(f).Fst1 = centerFreq(f)-sigma_f(f)-0.5;
  ConnectomeEnvelopeParams.sigBandpass(f).Fp1 = centerFreq(f)-sigma_f(f);
  ConnectomeEnvelopeParams.sigBandpass(f).Fp2 = centerFreq(f)+sigma_f(f);
  ConnectomeEnvelopeParams.sigBandpass(f).Fst2 = centerFreq(f)+sigma_f(f)+0.5;
  ConnectomeEnvelopeParams.sigBandpass(f).Ast1 = 40;
  ConnectomeEnvelopeParams.sigBandpass(f).Ap = 1;
  ConnectomeEnvelopeParams.sigBandpass(f).Ast2 = 40;
end

% ConnectomeEnvelopeParams.sigBandpass(1).Fst1 = 2.5;
% ConnectomeEnvelopeParams.sigBandpass(1).Fp1 = 3;
% ConnectomeEnvelopeParams.sigBandpass(1).Fp2 = 7;
% ConnectomeEnvelopeParams.sigBandpass(1).Fst2 = 7.5;
% ConnectomeEnvelopeParams.sigBandpass(1).Ast1 = 40;
% ConnectomeEnvelopeParams.sigBandpass(1).Ap = 1;
% ConnectomeEnvelopeParams.sigBandpass(1).Ast2 = 40;
% 
% ConnectomeEnvelopeParams.sigBandpass(2).Fst1 = 8.5;
% ConnectomeEnvelopeParams.sigBandpass(2).Fp1 = 9;
% ConnectomeEnvelopeParams.sigBandpass(2).Fp2 = 13;
% ConnectomeEnvelopeParams.sigBandpass(2).Fst2 = 13.5;
% ConnectomeEnvelopeParams.sigBandpass(2).Ast1 = 40;
% ConnectomeEnvelopeParams.sigBandpass(2).Ap = 1;
% ConnectomeEnvelopeParams.sigBandpass(2).Ast2 = 40;
% 
% ConnectomeEnvelopeParams.sigBandpass(3).Fst1 = 22.5;
% ConnectomeEnvelopeParams.sigBandpass(3).Fp1 = 23;
% ConnectomeEnvelopeParams.sigBandpass(3).Fp2 = 27;
% ConnectomeEnvelopeParams.sigBandpass(3).Fst2 = 27.5;
% ConnectomeEnvelopeParams.sigBandpass(3).Ast1 = 40;
% ConnectomeEnvelopeParams.sigBandpass(3).Ap = 1;
% ConnectomeEnvelopeParams.sigBandpass(3).Ast2 = 40;

ConnectomeEnvelopeParams.envLowpass.Fp = 0.5;
ConnectomeEnvelopeParams.envLowpass.Fst = 1;
ConnectomeEnvelopeParams.envLowpass.Ap = 1;
ConnectomeEnvelopeParams.envLowpass.Ast = 40;

ConnectomeEnvelopeParams.outFilenames = 'ConnectomeEnvelope';
ConnectomeEnvelopeParams.saveSourceRate = false;
ConnectomeEnvelopeParams.saveSourcePhase = false;
ConnectomeEnvelopeParams.saveSourceBP = false;
ConnectomeEnvelopeParams.saveEnvSig = false;
ConnectomeEnvelopeParams.saveEnvLP = false;

ConnectomeEnvelopeParams.deleteSimResult = false;

params.ConnectomeSim.ConnectomeEnvelope(1) = ConnectomeEnvelopeParams;
params.ConnectomeSim.ConnectomeEnvelope(1).applyLeadField = false;
params.ConnectomeSim.ConnectomeEnvelope(1).deleteSimResult = false;
params.ConnectomeSim.ConnectomeEnvelope(1).outFilenames = 'ConnectomeEnvelope';

params.ConnectomeSim.ConnectomeEnvelope(2) = ConnectomeEnvelopeParams;
params.ConnectomeSim.ConnectomeEnvelope(2).applyLeadField = true;
params.ConnectomeSim.ConnectomeEnvelope(2).applyLeadFieldMethod = 'perROI';
params.ConnectomeSim.ConnectomeEnvelope(2).deleteSimResult = false;
params.ConnectomeSim.ConnectomeEnvelope(2).outFilenames = 'ConnectomeEnvelopeLfPerROI';

params.ConnectomeSim.ConnectomeEnvelope(3) = ConnectomeEnvelopeParams;
params.ConnectomeSim.ConnectomeEnvelope(3).applyLeadField = true;
params.ConnectomeSim.ConnectomeEnvelope(3).applyLeadFieldMethod = 'perVertex';
params.ConnectomeSim.ConnectomeEnvelope(3).deleteSimResult = false;
params.ConnectomeSim.ConnectomeEnvelope(3).outFilenames = 'ConnectomeEnvelopeLfPerVertex';

paramsAll{1} = params;

[variableParams, paramComb] = getVariableParams(paramsAll{1},false);





clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG_ConnFC';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = 20;

params.ConnectomeEnvelopeReduce.onlyFCsim = false;
params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_ConnFC';
params.ConnectomeEnvelopeReduce.onlyCollectConnFC = true;
params.ConnectomeEnvelopeReduce.useEnvFreqAsParamVar = true;

paramsAll{2} = params;

clear eeg;
eeg.subj.Ids = [1:4 6:10 11:13 15 17:20];
eeg.subj.Avg = true;
eeg.day.Ids = 1;
eeg.day.Avg = false;
eeg.cond.Ids = 5:6;
eeg.cond.Avg = true;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;
params.Gridjob.runOnlyJobIds = [];

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {30, 25 ,27,28,29};
params.ConnectomeEnvelopeReduce.compareSC = false;
params.ConnectomeEnvelopeReduce.onlyFCsim = false;
params.ConnectomeEnvelopeReduce.useEnvFreqAsParamVar = true;

params.ConnectomeEnvelopeReduce.eeg = eeg;

params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;
params.ConnectomeEnvelopeReduce.sim.homotopic.Ids = [];
params.ConnectomeEnvelopeReduce.sim.homotopic.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG';
params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder = [];
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
paramsAll{3} = params;

%% lead field
clear params;
params.Gridjob.runLocal = false;
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

params.ConnectomeEnvelopeReduce.onlyFCsim = false;
params.ConnectomeEnvelopeReduce.outDirectory = {'CompareWithEEG_ConnFCLfPerROI','CompareWithEEG_ConnFCLfPerVertex'};
params.ConnectomeEnvelopeReduce.onlyCollectConnFC = true;
params.ConnectomeEnvelopeReduce.useEnvFreqAsParamVar = true;

paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = false;
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

params.ConnectomeEnvelopeReduce.sim.delay.Ids = [];%2;
params.ConnectomeEnvelopeReduce.sim.delay.Avg = false;
params.ConnectomeEnvelopeReduce.sim.v.Ids = [];
params.ConnectomeEnvelopeReduce.sim.v.Avg = false;

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
paramsAll{5} = params;


clear params;
gridjobs = Gridjob(paramsAll(3));
start(gridjobs);

