clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeSim';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 2; % 0=datasimu from Arnaud, 1=Bastian, 2=BastianNew
params.ConnectomeSim.subjId = -1;

params.ConnectomeSim.homotopic = 0;
params.ConnectomeSim.normRow = 1;
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
params.ConnectomeSim.k=700;
params.ConnectomeSim.v=1.75;
params.ConnectomeSim.delay=2; %typically 0.3-0.5 ms up to 2 ms
params.ConnectomeSim.t_max={620,920,1220};
params.ConnectomeSim.dt=0.0001;
params.ConnectomeSim.sampling=10;
params.ConnectomeSim.sig_n=0;
params.ConnectomeSim.d=0;
params.ConnectomeSim.verbose=true;
params.ConnectomeSim.statsRemoveInitialT = 20;
params.ConnectomeSim.outFilenames = 'ConnectomeSim';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeEnvelope';
params.Gridjob.combParallel = true;
params.Gridjob.initRandStreamWithSeed = 1;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.ConnectomeEnvelope.inFileRates = cellfun(@(x) ['ConnectomeSim/' num2str(x) 'SimResult.mat'], num2cell(1:3),'UniformOutput',false);
params.ConnectomeEnvelope.source_t_start = -Inf;
params.ConnectomeEnvelope.source_t_end = Inf;
params.ConnectomeEnvelope.saveSamples_t_start = [20 55];
params.ConnectomeEnvelope.saveSamples_t_end = [25 60];
params.ConnectomeEnvelope.env_t_start = 20;
params.ConnectomeEnvelope.env_t_end = {620,920,1220};
params.ConnectomeEnvelope.filtermethod = 'butter'; %or equiripple

params.ConnectomeEnvelope.sigBandpass(1).Fst1 = 2.5;
params.ConnectomeEnvelope.sigBandpass(1).Fp1 = 3;
params.ConnectomeEnvelope.sigBandpass(1).Fp2 = 7;
params.ConnectomeEnvelope.sigBandpass(1).Fst2 = 7.5;
params.ConnectomeEnvelope.sigBandpass(1).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(1).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(1).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(2).Fst1 = 8.5;
params.ConnectomeEnvelope.sigBandpass(2).Fp1 = 9;
params.ConnectomeEnvelope.sigBandpass(2).Fp2 = 13;
params.ConnectomeEnvelope.sigBandpass(2).Fst2 = 13.5;
params.ConnectomeEnvelope.sigBandpass(2).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(2).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(2).Ast2 = 40;

params.ConnectomeEnvelope.sigBandpass(3).Fst1 = 22.5;
params.ConnectomeEnvelope.sigBandpass(3).Fp1 = 23;
params.ConnectomeEnvelope.sigBandpass(3).Fp2 = 27;
params.ConnectomeEnvelope.sigBandpass(3).Fst2 = 27.5;
params.ConnectomeEnvelope.sigBandpass(3).Ast1 = 40;
params.ConnectomeEnvelope.sigBandpass(3).Ap = 1;
params.ConnectomeEnvelope.sigBandpass(3).Ast2 = 40;

params.ConnectomeEnvelope.envLowpass.Fp = 0.5;
params.ConnectomeEnvelope.envLowpass.Fst = 1;
params.ConnectomeEnvelope.envLowpass.Ap = 1;
params.ConnectomeEnvelope.envLowpass.Ast = 40;

params.ConnectomeEnvelope.outFilenames = 'ConnectomeEnvelope';
params.ConnectomeEnvelope.saveSourceRate = false;
params.ConnectomeEnvelope.saveSourcePhase = false;
params.ConnectomeEnvelope.saveSourceBP = false;
params.ConnectomeEnvelope.saveEnvSig = false;
params.ConnectomeEnvelope.saveEnvLP = false;
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeEnvelopeReduce';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegSubjIds = [1:4 7:10];
params.ConnectomeEnvelopeReduce.outDirectory = 'ConnectomeEnvelopeReduce';
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 14000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ConnectomeEnvelopeReduceDb2';
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '6';
params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'ConnectomeEnvelope';
params.ConnectomeEnvelopeReduce.eegDatabase = 2;
params.ConnectomeEnvelopeReduce.eegSubjIds = [1:4 7:10];
params.ConnectomeEnvelopeReduce.outDirectory = 'ConnectomeEnvelopeReduceDb2';
paramsAll{4} = params;

clear params;
gridjobs = Gridjob(paramsAll(4));
start(gridjobs);



