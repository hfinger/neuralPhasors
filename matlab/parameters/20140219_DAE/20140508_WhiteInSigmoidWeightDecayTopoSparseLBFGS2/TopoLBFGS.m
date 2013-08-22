clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'TopoLBFGS';
params.Gridjob.initRandStreamWithSeed = 12345;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '2';
params.Gridjob.combParallel = true;

params.Autoencoder.inActFolder = '../../20131217_ReLuDAE/20131220_MoreImages/labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'TopoLBFGS';
params.Autoencoder.continue = true;
params.Autoencoder.continueBatch = false;
params.Autoencoder.continueBatchInWeightsFolder = [];
params.Autoencoder.continueBatchInBackConnFilenames = [];
params.Autoencoder.continueBatchInForwConnFilenames = [];

% INIT PARAMS:
params.Autoencoder.inSamplesDims = [60 60 3 500]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 6;
params.Autoencoder.patchDimBackward = 6;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.useTiledConv = false;
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

params.Autoencoder.useRectifiedLinear = false;
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

params.Autoencoder.maskingNoiseFractionHidden =  0.5;        % fraction of inputs to set to 0
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
params.Autoencoder.gamma = 3;        % weight of topographic sparsity penalty term
params.Autoencoder.gammaPvalue = 0.5;
params.Autoencoder.topoNghFcn = {@(x) exp(-((1:x)-x/2).^2/4),@(x) exp(-((1:x)-x/2).^2/8)};
params.Autoencoder.topoNgh = {[9 9],[9 9]};
params.Autoencoder.topoGridDims = [10 10];
params.Autoencoder.topoPeriodicBoundary = [true true];
params.Autoencoder.topoEpsilon = 1e-2;
params.Autoencoder.mu = 0;        % weight of modular topographic sparsity penalty term
params.Autoencoder.muPvalue = 0.5;        % Psi_mu = mu * ( muFeatureKernel * (muConvKernel * y2.^2) + muEpsilon ) ^ muPvalue
params.Autoencoder.muConvKernel = 1; %fspecial('gaussian', 7, 2);
params.Autoencoder.muFeatureKernel = []; % if empty then no combination in feature dimension, otherwise it should be a matrix of size (hiddenSize,hiddenSize)
params.Autoencoder.muEpsilon = 1e-2;
params.Autoencoder.lambda = 1;%0.1;    % weight of L2-decay parameter
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
params.Autoencoder.saveintervalLog = 20;
params.Autoencoder.saveintervalpng = 1;
params.Autoencoder.saveintervalpngKeep = 100;
params.Autoencoder.savepng = true;
params.Autoencoder.reloadSaveinterval = 1;
params.Autoencoder.debugOutput = true;

% VALIDATION:
params.Autoencoder.validationSetsize = 20;
params.Autoencoder.validationInterval = Inf;
params.Autoencoder.validationIntervalLog = 1;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.validationSetImageIds = 1:10:1001;
params.Autoencoder.validationSkipGrad = false;

% ITERATION META-PARAMS:
params.Autoencoder.resetRandstreamEachReload = false;
params.Autoencoder.resetRandstreamEachIter = false;
params.Autoencoder.resetRandstreamEachEval = true;
params.Autoencoder.fixedBatches = false;
params.Autoencoder.batchsize = []; %[] means full batch
params.Autoencoder.numberOfPatchReloads =  200;
params.Autoencoder.numberOfImagesPerPatchReload = 200;
params.Autoencoder.useMinFuncGrad = false;

% params.minFuncGrad.Method = 'momentum';
% params.minFuncGrad.maxIter = 100;
% params.minFuncGrad.learnrate = 1e-3;
% params.minFuncGrad.momentum = 0.9;
% params.minFuncGrad.display = 'on';
% params.minFuncGrad.displayEvery = 1;
params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 20;
params.minFunc.display = 'on';
params.minFunc.optTol = 1e-6;
params.minFunc.progTol = 1e-14;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


