
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'AE100layer1PhaseTest';
params.PhaseSimulation.inActFolder = 'AE100layer1Act';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'AE100layer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'AE100layer1PhaseTest';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode1';
params.PhaseSimulation.weightAll = 30;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;

gridjobs = Gridjob(params);
start(gridjobs);
