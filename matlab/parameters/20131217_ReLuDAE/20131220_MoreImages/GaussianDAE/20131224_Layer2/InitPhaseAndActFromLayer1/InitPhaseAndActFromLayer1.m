clear paramsAll;

edges = 0:0.05:1;
edgesScaled = 0:0.1:10;

%% whitening without phase weighting
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'layer2InitActAfterWhiteningWithoutPhaseWeighting';
params.Gridjob.wc_host = '';
params.Gridjob.fhandleFinish = @() activityHist( 'layer2InitActAfterWhiteningWithoutPhaseWeighting', 'act1.mat', 'layer2InitActAfterWhiteningWithoutPhaseWeighting/stats.mat', edges, 100, true, edgesScaled );
params.ApplyWeights.inActFolder = '../NoPhaseFromLayer1/layer1ActRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = num2cell(1:50);
params.ApplyWeights.outActFolder = 'layer2InitActAfterWhiteningWithoutPhaseWeighting';
params.ApplyWeights.weightsFile = '../NoPhaseFromLayer1/layer1ActRectifiedWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;

%% whitening with phase weighting
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'layer2InitActAfterWhitening';
params.Gridjob.wc_host = '';
params.Gridjob.fhandleFinish = @() activityHist( 'layer2InitActAfterWhitening', 'act1.mat', 'layer2InitActAfterWhitening/stats.mat', edges, 100, true, edgesScaled );
params.ApplyWeightsWithTargetPhase.inActFolder = '../NoPhaseFromLayer1/layer1ActRectified';
params.ApplyWeightsWithTargetPhase.inActFilenames = 'act1.mat';
params.ApplyWeightsWithTargetPhase.inPhaseFolder = '../../20130726_Paper/Autoencoder/img50PhaseManySamplesFDR0p05maxdx32SmallDt';
params.ApplyWeightsWithTargetPhase.inPhaseFilenames = 'phaseIter30.mat';
params.ApplyWeightsWithTargetPhase.inTargetPhaseFolder = '../InitPhaseFromLayer1AbsW/layer2InitPhaseAfterWhitening'; 
params.ApplyWeightsWithTargetPhase.inTargetPhaseFilenames = 'phaseIter30.mat'; 
params.ApplyWeightsWithTargetPhase.fileid = num2cell(1:50);
params.ApplyWeightsWithTargetPhase.outActFolder = 'layer2InitActAfterWhitening';
params.ApplyWeightsWithTargetPhase.weightsFile = '../NoPhaseFromLayer1/layer1ActRectifiedWhiteningWeights/weights.mat';
params.ApplyWeightsWithTargetPhase.convType = 'same';
paramsAll{2} = params;


%% AE without phase weighting
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'layer2ActNotRectifiedWithoutPhaseWeighting';
params.Gridjob.fhandleFinish = @() activityHist( 'layer2ActNotRectifiedWithoutPhaseWeighting', 'act1.mat', 'layer2ActNotRectifiedWithoutPhaseWeighting/stats.mat', edges, 100, true, edgesScaled );
params.ApplyWeights.inActFolder = 'layer2InitActAfterWhiteningWithoutPhaseWeighting';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = num2cell(1:50);
params.ApplyWeights.outActFolder = 'layer2ActNotRectifiedWithoutPhaseWeighting';
params.ApplyWeights.weightsFile = '../NoPhaseFromLayer1/layer2AE100weights/2/forwConn.mat';
params.ApplyWeights.convType = 'same';
params.ApplyWeights.b = zeros(1,1,100); %delete AE bias for phase propagation
params.ApplyWeights.actFcn = @(x) x; %delete AE activation function for phase propagation
paramsAll{3} = params;

%% AE with phase weighting
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'layer2ActNotRectifiedPhaseWeighted';
params.Gridjob.wc_host = '';
params.Gridjob.fhandleFinish = @() activityHist( 'layer2ActNotRectifiedPhaseWeighted', 'act1.mat', 'layer2ActNotRectifiedPhaseWeighted/stats.mat', edges, 100, true, edgesScaled );
params.ApplyWeightsWithTargetPhase.inActFolder = 'layer2InitActAfterWhitening';
params.ApplyWeightsWithTargetPhase.inActFilenames = 'act1.mat';
params.ApplyWeightsWithTargetPhase.inPhaseFolder = '../InitPhaseFromLayer1AbsW/layer2InitPhaseAfterWhitening';
params.ApplyWeightsWithTargetPhase.inPhaseFilenames = 'phaseIter30.mat';
params.ApplyWeightsWithTargetPhase.inTargetPhaseFolder = '../InitPhaseFromLayer1AbsW/layer2InitPhaseAfterAE'; 
params.ApplyWeightsWithTargetPhase.inTargetPhaseFilenames = 'phaseIter30.mat'; 
params.ApplyWeightsWithTargetPhase.fileid = num2cell(1:50);
params.ApplyWeightsWithTargetPhase.outActFolder = 'layer2ActNotRectifiedPhaseWeighted';
params.ApplyWeightsWithTargetPhase.weightsFile = '../NoPhaseFromLayer1/layer2AE100weights/2/forwConn.mat';
params.ApplyWeightsWithTargetPhase.convType = 'same';
params.ApplyWeightsWithTargetPhase.b = zeros(1,1,100); %delete AE bias for phase propagation
params.ApplyWeightsWithTargetPhase.actFcn = @(x) x; %delete AE activation function for phase propagation
paramsAll{4} = params;

%% apply bias and activation function, but not W:
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActNotRectified';
params.ApplyWeights.inActFolder = 'layer2ActNotRectifiedPhaseWeighted';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:50;
params.ApplyWeights.outActFolder = 'layer2ActNotRectified';
params.ApplyWeights.weightsFile = '../NoPhaseFromLayer1/layer2AE100weights/2/forwConn.mat';
params.ApplyWeights.convType = 'same';
params.ApplyWeights.W = false;
paramsAll{5} = params;

%% continue as before:
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActRectified';
params.ApplyWeights.inActFolder = 'layer2ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:50;
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
params.ProbabilisticConn.exclDoubleConns = false;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0;
params.ProbabilisticConn.exclCorrWithPValueAbove = [];
params.ProbabilisticConn.FDRalpha = []; %0.05;
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
gridjobs = Gridjob(paramsAll(8:11));
start(gridjobs);


