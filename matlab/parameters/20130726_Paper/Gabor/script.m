%% create Gabor filter weights
rfSize = 12;
numOrient = 8;
order = 2;
sigma = 1.5;
aspectratio = 2;
integral = 0.1;
numChan = 3;

W = zeros(rfSize,rfSize,numChan,1,1,numChan*numOrient);
Wtmp = zeros(rfSize,rfSize,1,1,1,numOrient);
orientations=(0:numOrient-1)*pi/numOrient;
for k=1:length(orientations);
  Wtmp(:,:,1,1,1,k) = genSpatialKernel( rfSize, order, sigma, aspectratio, integral, orientations(k) );
end
for k=1:numChan
  W(:,:,k,1,1,(k-1)*numOrient+1:k*numOrient) = Wtmp;
end
conn.W = 5 * cat(6,W,-W); % add opposite colors
% so W in total consists of surround -5, center 10.5 and surround -5
% conn.inputSubtract = zeros(1,1,3);
% conn.inputSubtract(1,1,:) = [0.4255, 0.4409, 0.4324]; %subtract mean
% conn.inputScaling = zeros(1,1,3);
% conn.inputScaling(1,1,:) = [1/0.2465, 1/0.2432, 1/0.2563]; %standardize input
conn.inputSubsampling = 2;

mkdir('/net/store/nbp/phasesim/workdir/20130726_Paper/Gabor/GaborWeights')
save('/net/store/nbp/phasesim/workdir/20130726_Paper/Gabor/GaborWeights/forwConn.mat','-struct','conn')

%%
clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeInput';
params.LoadLabelMe.catName = '05june05_static_street_boston';
params.LoadLabelMe.fileid = 1:185;
params.LoadLabelMe.outActFolder = 'labelMeInput';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeWhiteningWeights';
params.Gridjob.requiremf = 12000;
params.Whitening.inActFolder = 'labelMeInput';
params.Whitening.outWeightsFolder = 'labelMeWhiteningWeights';
params.Whitening.inNumChannels = 3;
params.Whitening.maxCorrLengthDim1 = 64;
params.Whitening.maxCorrLengthDim2 = 64;
params.Whitening.convKernelDim1 = 51;
params.Whitening.convKernelDim2 = 51;
params.Whitening.reduceToConv = true;
params.Whitening.epsilon = 0.1;
params.Whitening.numPatches = 10000;
params.Whitening.matchFilenames = 'act.*.mat';
params.Whitening.borderBuffer = 0;
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeWhite';
params.ApplyWeights.inActFolder = 'labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:50;
params.ApplyWeights.outActFolder = 'labelMeWhite';
params.ApplyWeights.weightsFile = 'labelMeWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1ActNotRectified';
params.ApplyWeights.inActFolder = 'labelMeWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:50;
params.ApplyWeights.outActFolder = 'layer1ActNotRectified';
params.ApplyWeights.weightsFile = conn;
params.ApplyWeights.actFcn = @(x) 1 ./ (1 + exp(-x));
params.ApplyWeights.convType = 'same';
paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1ActRectified';
params.ApplyWeights.inActFolder = 'layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{5} = params;


clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1CovManySamples';
params.Gridjob.requiremf = 9000;
params.Gridjob.wc_host = '';
params.FeatureCovariance.inActFolder = 'layer1ActRectified';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = 1:50;
params.FeatureCovariance.outCovFolder = 'layer1CovManySamples';
params.FeatureCovariance.maxCovLengthDim1 = 32;
params.FeatureCovariance.maxCovLengthDim2 = 32;
params.FeatureCovariance.numSamplesPerImage = 500;
params.FeatureCovariance.borderBuffer = 0;
params.FeatureCovariance.saveCorrPvalue = true;
paramsAll{6} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1ConnManySamples';
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.ProbabilisticConn.inCorrFile = 'layer1CovManySamples/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1ConnManySamples';
params.ProbabilisticConn.numExc = 200;
params.ProbabilisticConn.numInh = 200;
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = true;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0;
params.ProbabilisticConn.exclCorrWithPValueAbove = [];
params.ProbabilisticConn.FDRalpha = 0.05;
params.ProbabilisticConn.maxdx = 32;
params.ProbabilisticConn.maxdy = 32;
paramsAll{7} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1PhaseAllManySamples';
params.Gridjob.wc_host = '';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:50);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer1ConnManySamples/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1PhaseAllManySamples';
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
paramsAll{8} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1EvalManySamples';
params.Gridjob.combParallel = true;
params.Gridjob.wc_host = '';
params.SegmentationEval.inActFolder = 'layer1ActRectified';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'layer1PhaseAllManySamples';
params.SegmentationEval.outSegEvalFolder = 'layer1EvalManySamplesBugFixed';
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
paramsAll{9} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1EvalManySamplesOnlyBorderEval';
params.Gridjob.combParallel = true;
params.Gridjob.wc_host = '';
params.SegmentationEval.inActFolder = 'layer1ActRectified';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'layer1PhaseAllManySamples';
params.SegmentationEval.outSegEvalFolder = 'layer1EvalManySamplesOnlyBorderEval';
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
paramsAll{10} = params;

clear params;
gridjobs = Gridjob(paramsAll(10));
start(gridjobs);
