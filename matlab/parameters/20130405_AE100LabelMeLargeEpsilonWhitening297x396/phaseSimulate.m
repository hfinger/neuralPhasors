% clear params;
% params.Gridjob.runLocal = true;
% params.Gridjob.requiremf = 12000;
% params.Gridjob.jobname = 'lMLayer1';
% params.ApplyWeights.inActFolder = 'labelMeWhite';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.outActFolder = 'AE100subsActLayer1';
% params.ApplyWeights.weightsFile = 'AE100subsVarWeights/6/forwConn.mat';
% params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), min(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
% params.ApplyWeights.convType = 'same';
% paramsAll{1} = params;

% clear params;
% params.Gridjob.runLocal = true;
% params.Gridjob.jobname = 'lMLayer1Cov';
% params.Gridjob.requiremf = 13000;
% params.FeatureCovariance.inActFolder = 'AE100subsActLayer1';
% params.FeatureCovariance.inActFilenames = 'act.*.mat';
% params.FeatureCovariance.fileid = [];
% params.FeatureCovariance.outCovFolder = 'labelMeLayer1Cov';
% params.FeatureCovariance.inNumFeatures = 100;
% params.FeatureCovariance.maxCovLengthDim1 = 13;
% params.FeatureCovariance.maxCovLengthDim2 = 13;
% params.FeatureCovariance.numSamplesPerImage = 20;
% params.FeatureCovariance.borderBuffer = 0;
% paramsAll{1} = params;
% 
% clear params;
% params.Gridjob.runLocal = true;
% params.Gridjob.jobname = 'lMLayer1Conn';
% params.Gridjob.requiremf = 13000;
% params.ProbabilisticConn.inCorrFile = 'labelMeLayer1Cov/patchCorr.mat';
% params.ProbabilisticConn.outWeightsFolder = 'labelMeLayer1Conn';
% params.ProbabilisticConn.numExc = 200;
% params.ProbabilisticConn.numInh = 200;
% params.ProbabilisticConn.exclSelfConn = true;
% params.ProbabilisticConn.useDiscretesample = true;
% params.ProbabilisticConn.exclDoubleConns = true;
% params.ProbabilisticConn.probRadiusFcn = [];
% paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'lMLayer1Phase';
params.PhaseSimulation.inActFolder = 'AE100subsActLayer1';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'labelMeLayer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'labelMeLayer1Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode1';
params.PhaseSimulation.weightAll = 10;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
