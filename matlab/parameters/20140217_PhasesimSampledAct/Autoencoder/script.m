clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1CovManySamples';
params.Gridjob.requiremf = 13000;
params.FeatureCovariance.inActFolder = '../../20130726_Paper/Autoencoder/layer1ActNotRectified';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = [];
params.FeatureCovariance.outCovFolder = 'layer1CovManySamples';
params.FeatureCovariance.maxCovLengthDim1 = 32;
params.FeatureCovariance.maxCovLengthDim2 = 32;
params.FeatureCovariance.numSamplesPerImage = 500;
params.FeatureCovariance.borderBuffer = 0;
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1ConnManySamplesFDR0p05maxdx32';
params.Gridjob.requiremf = 13000;
params.ProbabilisticConn.inCorrFile = 'layer1CovManySamples/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1ConnManySamplesFDR0p05maxdx32';
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
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'img50PhaseManySamplesFDR0p05maxdx32SmallDt';
params.PhaseSimulation.inActFolder = '../../20130726_Paper/Autoencoder/layer1ActNotRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:2);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer1ConnManySamplesFDR0p05maxdx32/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'phases';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 10;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
params.PhaseSimulation.sampleActRepeated = true;
paramsAll{3} = params;

clear params;
gridjobs = Gridjob(paramsAll(3));
start(gridjobs);


