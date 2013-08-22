
clear params;
params.Gridjob.runLocal = true;
params.Gridjob.jobname = 'layer1Conn0p15only100Conns';
params.Gridjob.requiremf = 3000;
params.ProbabilisticConn.inCorrFile = '../layer1Cov/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1Conn0p15only100Conns';
params.ProbabilisticConn.numExc = 100;
params.ProbabilisticConn.numInh = 100;
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = true;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0.15;
paramsAll{1} = params;


clear params;
gridjobs = Gridjob(paramsAll{1});
start(gridjobs);
