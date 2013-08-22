clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'ReLuDAE100weightsSubs2OptimizeAvgAct';
params.Gridjob.initRandStreamWithSeed = 12345;

params.Autoencoder.inActFolder = '../20131220_MoreImages/labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'ReLuDAE100weightsSubs2OptimizeAvgAct';
params.Autoencoder.continue = false;
params.Autoencoder.debugOutput = true;

params.Autoencoder.continueBatch = false;
params.Autoencoder.continueBatchInWeightsFolder = [];
params.Autoencoder.continueBatchInBackConnFilenames = [];
params.Autoencoder.continueBatchInForwConnFilenames = [];

params.Autoencoder.inSamplesDims = [60 60 3 2000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 8;
params.Autoencoder.patchDimBackward = 8;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = 2;
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = false;
params.Autoencoder.useNormMean = true;
params.Autoencoder.useRectifiedLinear = true;
params.Autoencoder.saveinterval = 10;
params.Autoencoder.saveintervalpng = 10;
params.Autoencoder.savepng = true;
params.Autoencoder.sparsityParam = 0.01;
params.Autoencoder.alpha = 1;        % weight of Autoencoder
params.Autoencoder.lambda = 0;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
params.Autoencoder.beta = {0.03, 0.1, 0.3, 1};         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0.5;

params.Autoencoder.numberOfPatchReloads = 50;
params.Autoencoder.numberOfImagesPerPatchReload = 200;

params.Autoencoder.batchsize = []; %[] means full batch
params.Autoencoder.fixedBatches = false;
params.Autoencoder.validationSetsize = 1000;
params.Autoencoder.validationInterval = 10;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.validationSetImageIds = [];

params.Autoencoder.resetRandstreamEachReload = false;
params.Autoencoder.resetRandstreamEachIter = false;
params.Autoencoder.resetRandstreamEachEval = true;

params.Autoencoder.useMinFuncGrad = false;

params.Autoencoder.topoNghFcn = [];
params.Autoencoder.topoNgh = [1];
params.Autoencoder.topoGridDims = [100];
params.Autoencoder.topoPeriodicBoundary = [false];
params.Autoencoder.topoEpsilon = 1e-2;

params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 20;
params.minFunc.display = 'on';
params.minFunc.optTol = 1e-6;
params.minFunc.progTol = 1e-14;

paramsAll{1} = params;



clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'ReLuDAE100weightsSubs2OptimizeAvgActWithoutActivationFunction';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'ReLuDAE100weightsSubs2OptimizeAvgActWithoutActivationFunction';
params.ApplyWeights.weightsFile = cellfun(@(x) ['ReLuDAE100weightsSubs2OptimizeAvgAct/' num2str(x) '/forwConn.mat'],num2cell(1:4),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn = @(x) x; 
params.ApplyWeights.actFcn2 = @(x) x;
paramsAll{2} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'ReLuDAE100weightsSubs2OptimizeAvgActWithActivationFunction';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'ReLuDAE100weightsSubs2OptimizeAvgActWithActivationFunction';
params.ApplyWeights.weightsFile = cellfun(@(x) ['ReLuDAE100weightsSubs2OptimizeAvgAct/' num2str(x) '/forwConn.mat'],num2cell(1:4),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;

gridjobs = Gridjob(paramsAll(2:3));
start(gridjobs);


