clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1ConnVarExcInh';
params.Gridjob.requiremf = 3000;
params.Gridjob.wc_host = '';
params.ProbabilisticConn.inCorrFile = 'layer1CovManySamples/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'layer1ConnVarExcInh';
params.ProbabilisticConn.numExc = {25,50,100,200,400,800,1600};
params.ProbabilisticConn.numInh = {25,50,100,200,400,800,1600};
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = true;
params.ProbabilisticConn.probRadiusFcn = [];
params.ProbabilisticConn.exclCorrBelow = 0;
params.ProbabilisticConn.exclCorrWithPValueAbove = [];
params.ProbabilisticConn.FDRalpha = 0.05;
params.ProbabilisticConn.maxdx = 32;
params.ProbabilisticConn.maxdy = 32;
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1PhaseVarExcInh';
params.Gridjob.wc_host = '';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = strcat('layer1ConnVarExcInh/', cellfun(@num2str,num2cell(1:49)','UniformOutput',false), '/weights.mat')';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1PhaseVarExcInh';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 10;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
