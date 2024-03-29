clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActRectified';
params.ApplyWeights.inActFolder = '../../20130726_Paper/Autoencoder/layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1ActRectifiedWhiteningWeights';
params.Gridjob.requiremf = 12000;
params.Gridjob.wc_host = '';
params.Whitening.inActFolder = 'layer1ActRectified';
params.Whitening.outWeightsFolder = 'layer1ActRectifiedWhiteningWeights';
params.Whitening.inNumChannels = 100;
params.Whitening.maxCorrLengthDim1 = 10;
params.Whitening.maxCorrLengthDim2 = 10;
params.Whitening.convKernelDim1 = 7;
params.Whitening.convKernelDim2 = 7;
params.Whitening.reduceToConv = true;
params.Whitening.epsilon = 0.1;
params.Whitening.numPatches = 1000;
params.Whitening.matchFilenames = 'act.*.mat';
params.Whitening.borderBuffer = 0;
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'layer1ActRectifiedWhite';
params.Gridjob.wc_host = '';
params.ApplyWeights.inActFolder = 'layer1ActRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = num2cell(1:185);
params.ApplyWeights.outActFolder = 'layer1ActRectifiedWhite';
params.ApplyWeights.weightsFile = 'layer1ActRectifiedWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 9000;
params.Gridjob.jobname = 'layer2AE100weights';
params.Gridjob.wc_host = '';
params.Autoencoder.inActFolder = 'layer1ActRectifiedWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'layer2AE100weights';
params.Autoencoder.continue = false;
params.Autoencoder.inSamplesDims = {[30 30 100 2000],[30 30 100 4000]}; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 6;
params.Autoencoder.patchDimBackward = 6;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = [];
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = false;
params.Autoencoder.useRectifiedLinear = false;
params.Autoencoder.saveinterval = 5;
params.Autoencoder.saveintervalpng = 5;
params.Autoencoder.savepng = false;
params.Autoencoder.sparsityParam = 0.035;
params.Autoencoder.alpha = 0.01/3;        % weight of Autoencoder
params.Autoencoder.lambda = 0.001;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
params.Autoencoder.beta = 0.3;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0;
params.Autoencoder.batchsize = 1000; %[] means full batch
params.Autoencoder.fixedBatches = false;
params.Autoencoder.validationSetsize = 1000;
params.Autoencoder.validationInterval = 10;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.resetRandstreamEachIter = true;
params.Autoencoder.useMinFuncGrad = false;
params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 500;
params.minFunc.display = 'on';
params.minFunc.optTol = 1e-6;
params.minFunc.progTol = 1e-14;
paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActNotRectified';
params.ApplyWeights.inActFolder = 'layer1ActRectifiedWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer2ActNotRectified';
params.ApplyWeights.weightsFile = 'layer2AE100weights/2/forwConn.mat';
params.ApplyWeights.convType = 'same';
paramsAll{5} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActRectified';
params.ApplyWeights.inActFolder = 'layer2ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer2ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{6} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer2Cov';
params.Gridjob.requiremf = 13000;
params.FeatureCovariance.inActFolder = 'layer2ActRectified';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = [];
params.FeatureCovariance.outCovFolder = 'layer2Cov';
params.FeatureCovariance.maxCovLengthDim1 = 32;
params.FeatureCovariance.maxCovLengthDim2 = 32;
params.FeatureCovariance.numSamplesPerImage = 50;
params.FeatureCovariance.borderBuffer = 0;
paramsAll{7} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer2Conn';
params.Gridjob.requiremf = 13000;
params.ProbabilisticConn.inCorrFile = 'layer2Cov/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer2Conn';
params.ProbabilisticConn.numExc = 200;
params.ProbabilisticConn.numInh = 200;
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = true;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0;
params.ProbabilisticConn.exclCorrWithPValueAbove = [];
params.ProbabilisticConn.FDRalpha = 0.05;
params.ProbabilisticConn.maxdx = 16;
params.ProbabilisticConn.maxdy = 16;
paramsAll{8} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'layer2Phase';
params.PhaseSimulation.inActFolder = 'layer2ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:50);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer2Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer2Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 3;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
paramsAll{9} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'segEval';
params.Gridjob.combParallel = true;
params.SegmentationEval.inActFolder = 'layer2ActRectified';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'layer2Phase';
params.SegmentationEval.outSegEvalFolder = 'segEval';
params.SegmentationEval.catName = '05june05_static_street_boston';
params.SegmentationEval.fileid = [1:50];
params.SegmentationEval.segSizeBinEdges = [exp(0:11)];
params.SegmentationEval.minPixelPerSeg = 1;
params.SegmentationEval.maxPixelPerSeg = Inf;
params.SegmentationEval.borderSize=0;
params.SegmentationEval.numPairs=10000;
params.SegmentationEval.numSubpop=100;
params.SegmentationEval.numSamplesInSubpop=1000;
params.SegmentationEval.phasePerPixel=false;
params.SegmentationEval.time={0,5,10,15,20,25,30};
params.SegmentationEval.initRandstream = true;
params.SegmentationEval.pairOnlyWithOneNonmatching = true;
params.SegmentationEval.verboseLevel = 0;
params.SegmentationEval.simOtherVarParam = [];
params.SegmentationEval.skipMissingImg = false;
params.SegmentationEval.borderEvalNumPerImg = 100;
params.SegmentationEval.borderEvalNghSize = 10;
params.SegmentationEval.borderAngleSigma = 3;
params.SegmentationEval.borderAngleUseActivityWeighting = false;
params.SegmentationEval.borderSyncUseActivityWeighting = true;
params.SegmentationEval.useHueInsteadPhase = false;
params.SegmentationEval.resizeToX = 200;
params.SegmentationEval.resizeToY = 150;
params.SegmentationEval.cropX = [];
params.SegmentationEval.cropY = [];
params.SegmentationEval.doCalcSegSync = true;
params.SegmentationEval.doCalcBorderSync = true;
params.SegmentationEval.doCalcBorderAngle = true;
paramsAll{10} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'segEvalOnlyBorder';
params.Gridjob.combParallel = true;
params.Gridjob.wc_host = '';
params.SegmentationEval.inActFolder = 'layer2ActRectified';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'layer2Phase';
params.SegmentationEval.outSegEvalFolder = 'segEvalOnlyBorder';
params.SegmentationEval.catName = '05june05_static_street_boston';
params.SegmentationEval.fileid = [1:50];
params.SegmentationEval.segSizeBinEdges = [exp(0:11)];
params.SegmentationEval.minPixelPerSeg = 36;
params.SegmentationEval.maxPixelPerSeg = 15000;
params.SegmentationEval.borderSize=0;
params.SegmentationEval.numPairs=10000;
params.SegmentationEval.numSubpop=100;
params.SegmentationEval.numSamplesInSubpop=1000;
params.SegmentationEval.phasePerPixel=false;
params.SegmentationEval.time={0,5,10,15,20,25,30};
params.SegmentationEval.initRandstream = true;
params.SegmentationEval.pairOnlyWithOneNonmatching = true;
params.SegmentationEval.verboseLevel = 0;
params.SegmentationEval.simOtherVarParam = [];
params.SegmentationEval.skipMissingImg = false;
params.SegmentationEval.borderEvalNumPerImg = 100;
params.SegmentationEval.borderEvalNghSize = 10;
params.SegmentationEval.borderAngleSigma = 3;
params.SegmentationEval.borderAngleUseActivityWeighting = false;
params.SegmentationEval.borderSyncUseActivityWeighting = true;
params.SegmentationEval.useHueInsteadPhase = false;
params.SegmentationEval.resizeToX = 200;
params.SegmentationEval.resizeToY = 150;
params.SegmentationEval.cropX = [];
params.SegmentationEval.cropY = [];
params.SegmentationEval.doCalcSegSync = false;
params.SegmentationEval.doCalcBorderSync = true;
params.SegmentationEval.doCalcBorderAngle = true;
paramsAll{11} = params;

clear params;
gridjobs = Gridjob(paramsAll(11));
start(gridjobs);


