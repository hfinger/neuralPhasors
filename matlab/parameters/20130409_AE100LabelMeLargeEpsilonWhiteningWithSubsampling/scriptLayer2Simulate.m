clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActNotRectified';
params.ApplyWeights.inActFolder = 'layer1ActRectifiedWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer2ActNotRectified';
params.ApplyWeights.weightsFile = 'temp_layer2AE100weights/forwConnIter80.mat';
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
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), min(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
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
paramsAll{8} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'layer2Phase';
params.PhaseSimulation.inActFolder = 'layer2ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
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
params.PhaseSimulation.weightAll = {10,30,100};
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;
paramsAll{9} = params;

clear params;
gridjobs = Gridjob(paramsAll(5:9));
start(gridjobs);


