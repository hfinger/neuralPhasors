clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ReLuDAE';
params.Gridjob.initRandStreamWithSeed = 1;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.Autoencoder.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'ReLuDAE';
params.Autoencoder.continue = true;
params.Autoencoder.continueBatch = false;
params.Autoencoder.continueBatchInWeightsFolder = [];
params.Autoencoder.continueBatchInBackConnFilenames = [];
params.Autoencoder.continueBatchInForwConnFilenames = [];

% INIT PARAMS:
params.Autoencoder.inSamplesDims = [120 120 3 2100]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 3;
params.Autoencoder.patchDimBackward = 3;
params.Autoencoder.hiddenSize = 25;
params.Autoencoder.useTiledConv = true;
params.Autoencoder.inputSubsampling = 5;
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
params.Autoencoder.nu = {0, 1, 10};
params.Autoencoder.nuEMAconst = 0.95;

nuCovWeightMat = abs(bsxfun( @minus, 1:25, (1:25)'));
nuCovWeightMat = reshape(nuCovWeightMat,[1 1 25 1 1 25]);
nuCovWeightMatFull = zeros([2 2 25 2 2 25]);
nuCovWeightMatFull(1,1,:,1,1,:) = nuCovWeightMat;
nuCovWeightMatFull(1,2,:,1,2,:) = nuCovWeightMat;
nuCovWeightMatFull(2,1,:,2,1,:) = nuCovWeightMat;
nuCovWeightMatFull(2,2,:,2,2,:) = nuCovWeightMat;

params.Autoencoder.nuCovWeightMat = reshape(nuCovWeightMatFull,[100 100]);
% params.Autoencoder.nuCovWeightMat = diag(ones(1,99),1)+diag(ones(1,99),-1)+diag(ones(1,100),0);

params.Autoencoder.lambda = 0;    % weight of L2-decay parameter
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
params.Autoencoder.saveintervalBackConn = 500;
params.Autoencoder.saveintervalTrainInfo = Inf;
params.Autoencoder.saveintervalLog = 100;
params.Autoencoder.saveintervalpng = 100;
params.Autoencoder.saveintervalpngKeep = 100;
params.Autoencoder.savepng = true;
params.Autoencoder.plotPermute = []%@(W) reshape(W,[size(W,1) size(W,2) size(W,3) 2 2 25]);
params.Autoencoder.reloadSaveinterval = 1;
params.Autoencoder.debugOutput = true;

% VALIDATION:
params.Autoencoder.validationSetsize = 100;
params.Autoencoder.validationInterval = Inf;
params.Autoencoder.validationIntervalLog = 100;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.validationSetImageIds = 1:10:1001;
params.Autoencoder.validationSkipGrad = false;

% ITERATION META-PARAMS:
params.Autoencoder.resetRandstreamEachReload = false;
params.Autoencoder.resetRandstreamEachIter = false;
params.Autoencoder.resetRandstreamEachEval = false;
params.Autoencoder.fixedBatches = false;
params.Autoencoder.batchsize = 10; %[] means full batch
params.Autoencoder.numberOfPatchReloads = 5;
params.Autoencoder.numberOfImagesPerPatchReload = 400;
params.Autoencoder.useMinFuncGrad = true;
params.Autoencoder.initialItersWithoutLearning = 20;

params.minFuncGrad.Method = 'momentum';
params.minFuncGrad.maxIter = 200;
params.minFuncGrad.learnrate = 1e-2;
params.minFuncGrad.momentum = 0.9;
params.minFuncGrad.display = 'on';
params.minFuncGrad.displayEvery = 10;

paramsAll{1} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActNotRectified';
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer1ActNotRectified';
params.ApplyWeights.weightsFile = 'ReLuDAE/forwConn.mat';
params.ApplyWeights.convType = 'same';
params.ApplyWeights.shiftOutputdims = true;
paramsAll{2} = params;


clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1CovManySamples';
params.Gridjob.wc_host = [];
params.Gridjob.requiremf = 13000;
params.Gridjob.requiredThreads = '4';
params.FeatureCovariance.inActFolder = 'layer1ActNotRectified';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = [];
params.FeatureCovariance.outCovFolder = 'layer1CovManySamples';
params.FeatureCovariance.maxCovLengthDim1 = 8;
params.FeatureCovariance.maxCovLengthDim2 = 8;
params.FeatureCovariance.numSamplesPerImage = 50;
params.FeatureCovariance.borderBuffer = 0;
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1ConnManySamplesFDR0p05maxdx8';
params.Gridjob.requiremf = 13000;
params.Gridjob.wc_host = [];
params.ProbabilisticConn.inCorrFile = 'layer1CovManySamples/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1ConnManySamplesFDR0p05maxdx8';
params.ProbabilisticConn.numExc = 200;
params.ProbabilisticConn.numInh = 200;
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = false;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0;
params.ProbabilisticConn.exclCorrWithPValueAbove = [];
params.ProbabilisticConn.FDRalpha = 0.05;
params.ProbabilisticConn.maxdx = 8;
params.ProbabilisticConn.maxdy = 8;
paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'phase';
params.Gridjob.wc_host = [];
params.PhaseSimulation.inActFolder = 'layer1ActNotRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer1ConnManySamplesFDR0p05maxdx8/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = {0.03, 0.1, 0.3};
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 2;
params.PhaseSimulation.saveintervalMeanPhase = 6;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 2;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
paramsAll{5} = params;

clear params;
gridjobs = Gridjob(paramsAll{1});
start(gridjobs);


