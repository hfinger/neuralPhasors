clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1CovManySamples2';
params.Gridjob.requiremf = 13000;
params.FeatureCovariance.inActFolder = 'layer1ActRectified';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = [];
params.FeatureCovariance.outCovFolder = 'layer1CovManySamples2';
params.FeatureCovariance.maxCovLengthDim1 = 32;
params.FeatureCovariance.maxCovLengthDim2 = 32;
params.FeatureCovariance.numSamplesPerImage = 1000;
params.FeatureCovariance.borderBuffer = 0;
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.jobname = 'layer1ConnManySamples2FDR0p05maxdx32';
params.Gridjob.requiremf = 13000;
params.ProbabilisticConn.inCorrFile = 'layer1CovManySamples2/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1ConnManySamples2FDR0p05maxdx32';
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
paramsAll{4} = params;


clear params;
gridjobs = Gridjob(paramsAll(3:4));
start(gridjobs);
