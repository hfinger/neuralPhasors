%% create Gabor filter weights
rfSize = 12;
numOrient = 8;
order = 2;
sigma = 1.5;
aspectratio = 2;
integral = 0.02;
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
% so W in total consists of surround -5, center 10.1 and surround -5
% conn.inputSubtract = zeros(1,1,3);
% conn.inputSubtract(1,1,:) = [0.4255, 0.4409, 0.4324]; %subtract mean
% conn.inputScaling = zeros(1,1,3);
% conn.inputScaling(1,1,:) = [1/0.2465, 1/0.2432, 1/0.2563]; %standardize input
conn.inputSubsampling = 2;

%%
clear paramsAll;

% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.requiremf = 3000;
% params.Gridjob.jobname = 'layer1ActNotRectified';
% params.ApplyWeights.inActFolder = '../20130409_AE100LabelMeLargeEpsilonWhiteningWithSubsampling/labelMeWhite';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.outActFolder = 'layer1ActNotRectified';
% params.ApplyWeights.weightsFile = conn;
% params.ApplyWeights.actFcn = @(x) 1 ./ (1 + exp(-x));
% params.ApplyWeights.convType = 'same';
% paramsAll{1} = params;
% 
% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.requiremf = 3000;
% params.Gridjob.jobname = 'layer1ActRectified';
% params.ApplyWeights.inActFolder = 'layer1ActNotRectified';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.outActFolder = 'layer1ActRectified';
% params.ApplyWeights.weightsFile = [];
% params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), min(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
% params.ApplyWeights.convType = 'same';
% paramsAll{2} = params;
% 
% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.jobname = 'layer1Cov';
% params.Gridjob.requiremf = 3000;
% params.FeatureCovariance.inActFolder = 'layer1ActRectified';
% params.FeatureCovariance.inActFilenames = 'act.*.mat';
% params.FeatureCovariance.fileid = 1:50;
% params.FeatureCovariance.outCovFolder = 'layer1Cov';
% params.FeatureCovariance.maxCovLengthDim1 = 18;
% params.FeatureCovariance.maxCovLengthDim2 = 18;
% params.FeatureCovariance.numSamplesPerImage = 100;
% params.FeatureCovariance.borderBuffer = 0;
% params.FeatureCovariance.saveCorrPvalue = true;
% paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1ConnFDR0p05';
params.Gridjob.requiremf = 3000;
params.ProbabilisticConn.inCorrFile = 'layer1Cov/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1ConnFDR0p05';
params.ProbabilisticConn.numExc = 200;
params.ProbabilisticConn.numInh = 200;
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = false;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0;
params.ProbabilisticConn.exclCorrWithPValueAbove = [];
params.ProbabilisticConn.FDRalpha = 0.05;
params.ProbabilisticConn.maxdx = 18;
params.ProbabilisticConn.maxdy = 18;
paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1PhaseFDR0p05';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer1ConnFDR0p05/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1PhaseFDR0p05';
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
params.PhaseSimulation.maxdphase = 0.5;
paramsAll{5} = params;

clear params;
gridjobs = Gridjob(paramsAll(4:5));
start(gridjobs);
