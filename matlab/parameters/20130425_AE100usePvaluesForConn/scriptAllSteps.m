clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1Cov';
params.Gridjob.requiremf = 13000;
params.FeatureCovariance.inActFolder = '../20130409_AE100LabelMeLargeEpsilonWhiteningWithSubsampling/layer1ActRectified';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = 1:50;
params.FeatureCovariance.outCovFolder = 'layer1Cov';
params.FeatureCovariance.maxCovLengthDim1 = 18;
params.FeatureCovariance.maxCovLengthDim2 = 18;
params.FeatureCovariance.numSamplesPerImage = 100;
params.FeatureCovariance.borderBuffer = 0;
params.FeatureCovariance.saveCorrPvalue = true;
paramsAll{1} = params;

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
params.ProbabilisticConn.exclDoubleConns = false;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0;
params.ProbabilisticConn.exclCorrWithPValueAbove = 0.05 / (100*(18*2+1)^2);
params.ProbabilisticConn.maxdx = 18;
params.ProbabilisticConn.maxdy = 18;
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'layer1Phase';
params.PhaseSimulation.inActFolder = '../20130409_AE100LabelMeLargeEpsilonWhiteningWithSubsampling/layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 50;
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
paramsAll{3} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


