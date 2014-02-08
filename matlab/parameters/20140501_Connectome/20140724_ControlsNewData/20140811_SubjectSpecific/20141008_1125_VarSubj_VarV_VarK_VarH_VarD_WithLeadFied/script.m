clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 6000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 2; % 0=datasimu from Arnaud, 1=Bastian, 2=BastianNew
params.ConnectomeSim.subjId = num2cell([1:4 6:10]);

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = num2cell(0:0.04:0.12);
params.ConnectomeSim.normRow = 1;
params.ConnectomeSim.shuffleSC = false;
params.ConnectomeSim.shuffleD = false;

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
params.ConnectomeSim.k=num2cell(100*2.^(0:1:3));
params.ConnectomeSim.v=num2cell(exp(-1:0.5:1));
params.ConnectomeSim.delay=num2cell(0:0.5:1); %typically 0.3-0.5 ms up to 2 ms
params.ConnectomeSim.t_max={290,290};
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 10;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';

ConnectomeEnvelopeParams.inFileRates = [];
ConnectomeEnvelopeParams.source_t_start = -Inf;
ConnectomeEnvelopeParams.source_t_end = Inf;
ConnectomeEnvelopeParams.saveSamples_t_start = [20 55];
ConnectomeEnvelopeParams.saveSamples_t_end = [25 60];
ConnectomeEnvelopeParams.env_t_start = 10;
ConnectomeEnvelopeParams.env_t_end = Inf;
ConnectomeEnvelopeParams.filtermethod = 'butter'; %or equiripple

ConnectomeEnvelopeParams.applyLeadField = true;
ConnectomeEnvelopeParams.applyLeadFieldMethod = 'maxpower';

ConnectomeEnvelopeParams.sigBandpass(1).Fst1 = 2.5;
ConnectomeEnvelopeParams.sigBandpass(1).Fp1 = 3;
ConnectomeEnvelopeParams.sigBandpass(1).Fp2 = 7;
ConnectomeEnvelopeParams.sigBandpass(1).Fst2 = 7.5;
ConnectomeEnvelopeParams.sigBandpass(1).Ast1 = 40;
ConnectomeEnvelopeParams.sigBandpass(1).Ap = 1;
ConnectomeEnvelopeParams.sigBandpass(1).Ast2 = 40;

ConnectomeEnvelopeParams.sigBandpass(2).Fst1 = 8.5;
ConnectomeEnvelopeParams.sigBandpass(2).Fp1 = 9;
ConnectomeEnvelopeParams.sigBandpass(2).Fp2 = 13;
ConnectomeEnvelopeParams.sigBandpass(2).Fst2 = 13.5;
ConnectomeEnvelopeParams.sigBandpass(2).Ast1 = 40;
ConnectomeEnvelopeParams.sigBandpass(2).Ap = 1;
ConnectomeEnvelopeParams.sigBandpass(2).Ast2 = 40;

ConnectomeEnvelopeParams.sigBandpass(3).Fst1 = 22.5;
ConnectomeEnvelopeParams.sigBandpass(3).Fp1 = 23;
ConnectomeEnvelopeParams.sigBandpass(3).Fp2 = 27;
ConnectomeEnvelopeParams.sigBandpass(3).Fst2 = 27.5;
ConnectomeEnvelopeParams.sigBandpass(3).Ast1 = 40;
ConnectomeEnvelopeParams.sigBandpass(3).Ap = 1;
ConnectomeEnvelopeParams.sigBandpass(3).Ast2 = 40;

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
params.ConnectomeSim.ConnectomeEnvelope(2).applyLeadFieldMethod = 'sum';
params.ConnectomeSim.ConnectomeEnvelope(2).deleteSimResult = false;
params.ConnectomeSim.ConnectomeEnvelope(2).outFilenames = 'ConnectomeEnvelopeLfSum';

params.ConnectomeSim.ConnectomeEnvelope(3) = ConnectomeEnvelopeParams;
params.ConnectomeSim.ConnectomeEnvelope(3).applyLeadField = true;
params.ConnectomeSim.ConnectomeEnvelope(3).applyLeadFieldMethod = 'eucl';
params.ConnectomeSim.ConnectomeEnvelope(3).deleteSimResult = false;
params.ConnectomeSim.ConnectomeEnvelope(3).outFilenames = 'ConnectomeEnvelopeLfEucl';

params.ConnectomeSim.ConnectomeEnvelope(4) = ConnectomeEnvelopeParams;
params.ConnectomeSim.ConnectomeEnvelope(4).applyLeadField = true;
params.ConnectomeSim.ConnectomeEnvelope(4).applyLeadFieldMethod = 'maxpower';
params.ConnectomeSim.ConnectomeEnvelope(4).deleteSimResult = true;
params.ConnectomeSim.ConnectomeEnvelope(4).outFilenames = 'ConnectomeEnvelopeLfMaxpower';


paramsAll{1} = params;

[variableParams, paramComb] = getVariableParams(paramsAll{1},false);


clear eeg;
eeg.subj.Ids = [1:4 6:10];
eeg.subj.Avg = false;
eeg.day.Ids = [];
eeg.day.Avg = true;
eeg.cond.Ids = 5:6;
eeg.cond.Avg = true;

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
params.ConnectomeEnvelopeReduce.eegDatabase = 10;

params.ConnectomeEnvelopeReduce.onlyFCsim = false;
params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_ConnFC';
params.ConnectomeEnvelopeReduce.onlyCollectConnFC = true;

paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = '!ramsauer.ikw.uni-osnabrueck.de';
params.Gridjob.jobname = 'CompareWithEEG_onlyMatchWithAicohcoh';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.Gridjob.combParallel = true;

params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = {10,11,12,13};
params.ConnectomeEnvelopeReduce.compareSC = false;
params.ConnectomeEnvelopeReduce.onlyFCsim = false;

params.ConnectomeEnvelopeReduce.eeg = eeg;

params.ConnectomeEnvelopeReduce.sim.t_max.Ids = [];
params.ConnectomeEnvelopeReduce.sim.t_max.Avg = true;
params.ConnectomeEnvelopeReduce.sim.delay.Ids = [];%2;
params.ConnectomeEnvelopeReduce.sim.delay.Avg = false;
params.ConnectomeEnvelopeReduce.sim.k.Ids = [];
params.ConnectomeEnvelopeReduce.sim.k.Avg = false;
params.ConnectomeEnvelopeReduce.sim.v.Ids = [];
params.ConnectomeEnvelopeReduce.sim.v.Avg = false;

params.ConnectomeEnvelopeReduce.outDirectory = 'CompareWithEEG_onlyMatchWithAicohcoh';
params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder = [];%cellfun(@(x) ['db' num2str(x)],params.ConnectomeEnvelopeReduce.eegDatabase,'UniformOutput',false);

params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'match';
params.ConnectomeEnvelopeReduce.results = struct();

params.ConnectomeEnvelopeReduce.calcCompareAicohWithCoh = true;

params.ConnectomeEnvelopeReduce.reloadConnFC = '../CompareWithEEG_ConnFC/ConnFC.mat';
params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
params.ConnectomeEnvelopeReduce.permutePlotDims = [2 3 1 4 5];
params.ConnectomeEnvelopeReduce.plotPermTests = false;
params.ConnectomeEnvelopeReduce.excludePlotFieldsRegexp = 'overFreq';
params.ConnectomeEnvelopeReduce.deletePlotFolder = false;
paramsAll{3} = params;


clear params;
gridjobs = Gridjob(paramsAll(1:3));
start(gridjobs);

