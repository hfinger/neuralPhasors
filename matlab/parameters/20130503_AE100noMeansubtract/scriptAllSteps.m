clear paramsAll;

% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.jobname = 'labelMeInput';
% params.LoadLabelMe.catName = '05june05_static_street_boston';
% params.LoadLabelMe.fileid = 1:185;
% params.LoadLabelMe.outActFolder = 'labelMeInput';
% paramsAll{1} = params;
% 
% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.jobname = 'labelMeWhiteningWeights';
% params.Gridjob.requiremf = 12000;
% params.Whitening.inActFolder = 'labelMeInput';
% params.Whitening.outWeightsFolder = 'labelMeWhiteningWeights';
% params.Whitening.inNumChannels = 3;
% params.Whitening.maxCorrLengthDim1 = 64;
% params.Whitening.maxCorrLengthDim2 = 64;
% params.Whitening.convKernelDim1 = 51;
% params.Whitening.convKernelDim2 = 51;
% params.Whitening.reduceToConv = true;
% params.Whitening.epsilon = 0.1;
% params.Whitening.numPatches = 10000;
% params.Whitening.matchFilenames = 'act.*.mat';
% params.Whitening.borderBuffer = 0;
% paramsAll{2} = params;
% 
% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.requiremf = 13000;
% params.Gridjob.jobname = 'labelMeWhite';
% params.ApplyWeights.inActFolder = 'labelMeInput';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.fileid = num2cell(1:185);
% params.ApplyWeights.outActFolder = 'labelMeWhite';
% params.ApplyWeights.weightsFile = 'labelMeWhiteningWeights/weights.mat';
% params.ApplyWeights.convType = 'same';
% paramsAll{3} = params;
% 
% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.requiremf = 5000;
% params.Gridjob.jobname = 'labelMeAE100weights';
% params.Autoencoder.inActFolder = 'labelMeWhite';
% params.Autoencoder.inActFilenames = 'act.*.mat';
% params.Autoencoder.outWeightsFolder = 'labelMeAE100weights';
% params.Autoencoder.continue = false;
% params.Autoencoder.inSamplesDims = [60 60 3 2000]; % [x,y,#features,#samples]
% params.Autoencoder.inSamplesBorderBuffer = 0;
% params.Autoencoder.patchDimForward = 6;
% params.Autoencoder.patchDimBackward = 6;
% params.Autoencoder.hiddenSize = 100;
% params.Autoencoder.inputSubsampling = 2;
% params.Autoencoder.forwInitMaxWeight = 0.01;
% params.Autoencoder.forwInitScaleDistFcn = [];
% params.Autoencoder.backInitMaxWeight = 0.01;
% params.Autoencoder.backInitScaleDistFcn = [];
% params.Autoencoder.useSoftmax = false;
% params.Autoencoder.useRectifiedLinear = false;
% params.Autoencoder.saveinterval = 10;
% params.Autoencoder.saveintervalpng = 5;
% params.Autoencoder.savepng = true;
% params.Autoencoder.sparsityParam = 0.035;
% params.Autoencoder.alpha = 0.01/3;        % weight of Autoencoder
% params.Autoencoder.lambda = 0.001;    % weight of L2-decay parameter
% params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
% params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
% params.Autoencoder.beta = 0.3;         % weight of sparsity penalty term
% params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
% params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
% params.Autoencoder.gaussianNoiseSigmaInput = 0;
% params.Autoencoder.batchsize = 1000; %[] means full batch
% params.Autoencoder.fixedBatches = false;
% params.Autoencoder.validationSetsize = 1000;
% params.Autoencoder.validationInterval = 10;
% params.Autoencoder.validationSetIds = [];
% params.Autoencoder.resetRandstreamEachIter = true;
% params.Autoencoder.useMinFuncGrad = false;
% params.minFunc.Method = 'lbfgs';
% params.minFunc.maxIter = 500;
% params.minFunc.display = 'on';
% params.minFunc.optTol = 1e-6;
% params.minFunc.progTol = 1e-14;
% paramsAll{4} = params;
% 
% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.requiremf = 12000;
% params.Gridjob.jobname = 'layer1ActNotRectified';
% params.ApplyWeights.inActFolder = 'labelMeWhite';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.fileid = 1:50;
% params.ApplyWeights.outActFolder = 'layer1ActNotRectified';
% params.ApplyWeights.weightsFile = 'labelMeAE100weights/forwConn.mat';
% params.ApplyWeights.convType = 'same';
% paramsAll{5} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActRectified';
params.ApplyWeights.inActFolder = '../20130409_AE100LabelMeLargeEpsilonWhiteningWithSubsampling/layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:50;
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) bsxfun(@rdivide,x,0.001+sum(x,3)); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{6} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1Cov';
params.Gridjob.requiremf = 13000;
params.FeatureCovariance.inActFolder = 'layer1ActRectified';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = [];
params.FeatureCovariance.outCovFolder = 'layer1Cov';
params.FeatureCovariance.maxCovLengthDim1 = 32;
params.FeatureCovariance.maxCovLengthDim2 = 32;
params.FeatureCovariance.numSamplesPerImage = 50;
params.FeatureCovariance.borderBuffer = 0;
paramsAll{7} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1Conn';
params.Gridjob.requiremf = 13000;
params.ProbabilisticConn.inCorrFile = 'layer1Cov/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1Conn';
params.ProbabilisticConn.numExc = 200;
params.ProbabilisticConn.numInh = 200;
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = true;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.FDRalpha = 0.01;
params.ProbabilisticConn.maxdx = 32;
params.ProbabilisticConn.maxdy = 32;
paramsAll{8} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'layer1Phase';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 30;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;
paramsAll{9} = params;

clear params;
gridjobs = Gridjob(paramsAll(6:9));
start(gridjobs);


