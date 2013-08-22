

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'gaborSimVideo';
params.Gridjob.wc_host = '';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'layer1ConnManySamples/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'gaborSimVideo';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 100;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 3;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
paramsAll{5} = params;


clear params;
gridjobs = Gridjob(paramsAll(5));
start(gridjobs);
