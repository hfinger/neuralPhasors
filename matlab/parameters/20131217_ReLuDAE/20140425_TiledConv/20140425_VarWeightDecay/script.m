clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ReLuDAEVarWD';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.Gridjob.remoteStart = true;

params.Autoencoder.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'ReLuDAEVarWD';
params.Autoencoder.continue = true;
params.Autoencoder.continueBatch = false;
params.Autoencoder.continueBatchInWeightsFolder = [];
params.Autoencoder.continueBatchInBackConnFilenames = [];
params.Autoencoder.continueBatchInForwConnFilenames = [];

% INIT PARAMS:
params.Autoencoder.inSamplesDims = [120 120 3 5500]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 7;
params.Autoencoder.patchDimBackward = 7;
params.Autoencoder.hiddenSize = 3;
params.Autoencoder.useTiledConv = true;
params.Autoencoder.inputSubsampling = 2;
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.initWGaussian = false;
params.Autoencoder.initWExponent = 1;

% PROCESSING STAGES:
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.maskingNoiseFractionRescale = true;        % if inputs are rescaled to the variance before masking
params.Autoencoder.gaussianNoiseSigmaInput = 0.75;
params.Autoencoder.gaussianNoiseCorrelated = false;
params.Autoencoder.gaussianNoiseCorrelatedUniform = false;

params.Autoencoder.useRectifiedLinear = true;
params.Autoencoder.useRectifiedLog = false; %instead of sigmoid
params.Autoencoder.useNoActFcn = false; %instead of sigmoid
params.Autoencoder.useBinaryLinearInterp = false;
params.Autoencoder.binaryLinearInterpScaling = 1;
params.Autoencoder.binaryLinearInterpEMAconst = 0.9;

params.Autoencoder.useSoftmax = false;
params.Autoencoder.useNormMean = false;
params.Autoencoder.normMeanEpsilon = 1e-6;

params.Autoencoder.hiddenLinearBinaryThresholdUnits = false; %these are not considered when calculating gradients
params.Autoencoder.hiddenSigmoidBinaryThresholdUnits = false; %these are not considered when calculating gradients

params.Autoencoder.maskingNoiseFractionHidden = 0;        % fraction of inputs to set to 0
params.Autoencoder.maskingNoiseFractionHiddenRescale = false;        % if inputs are rescaled to the variance before masking

params.Autoencoder.noBiasBack = false;

params.Autoencoder.reconstrSigmoid = false;
params.Autoencoder.reconstrRectifiedLinear = false;

% OPTIMIZATION PARAMS:
params.Autoencoder.alpha = 1;        % weight of Autoencoder
params.Autoencoder.useLinearAutoencError = false;
params.Autoencoder.useCrossEntropyError = false;
params.Autoencoder.beta = 0;         % weight of sparsity penalty term
params.Autoencoder.sparsityParam = 0.01;
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.gammaPvalue = 0.5;
params.Autoencoder.topoNghFcn = [];
params.Autoencoder.topoNgh = [1];
params.Autoencoder.topoGridDims = [100];
params.Autoencoder.topoPeriodicBoundary = [false];
params.Autoencoder.topoEpsilon = 1e-2;
params.Autoencoder.lambda = {0.03, 0.1, 0.3};    % weight of L2-decay parameter
params.Autoencoder.lambdaBackScale = 1;    % weight of L2-decay parameter scaling for backward conenctions
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.lambdaL1BackScale = 1;    % weight of L1-decay parameter
params.Autoencoder.fixL2WeightBack = false; %if true, remember to also set both backScale to 0
params.Autoencoder.fixL2WeightBackPoolPerHiddenNeuron = false;
params.Autoencoder.smoothLambdaL1WithEpsilon = 0 ; %if not zero then smooth it
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;

% LOGS AND SAVES:
params.Autoencoder.saveinterval = Inf;
params.Autoencoder.saveintervalForwConn = 500;
params.Autoencoder.saveintervalBackConn = Inf;
params.Autoencoder.saveintervalTrainInfo = Inf;
params.Autoencoder.saveintervalLog = 100;
params.Autoencoder.saveintervalpng = 100;
params.Autoencoder.savepng = true;
params.Autoencoder.reloadSaveinterval = 1;
params.Autoencoder.debugOutput = true;

% VALIDATION:
params.Autoencoder.validationSetsize = 500;
params.Autoencoder.validationInterval = Inf;
params.Autoencoder.validationIntervalLog = 100;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.validationSetImageIds = 1:10:1001;

% ITERATION META-PARAMS:
params.Autoencoder.resetRandstreamEachReload = false;
params.Autoencoder.resetRandstreamEachIter = false;
params.Autoencoder.resetRandstreamEachEval = false;
params.Autoencoder.fixedBatches = false;
params.Autoencoder.batchsize = 10; %[] means full batch
params.Autoencoder.numberOfPatchReloads = 100;
params.Autoencoder.numberOfImagesPerPatchReload = 500;
params.Autoencoder.useMinFuncGrad = true;

params.minFuncGrad.Method = 'momentum';
params.minFuncGrad.maxIter = 1000;
params.minFuncGrad.learnrate = 1e-2;
params.minFuncGrad.momentum = 0.9;
params.minFuncGrad.display = 'on';
params.minFuncGrad.displayEvery = 10;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


