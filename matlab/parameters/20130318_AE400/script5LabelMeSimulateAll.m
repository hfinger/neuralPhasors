clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'lMLayer1PhaseAll';
params.PhaseSimulation.inActFolder = 'labelMeActLayer1';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:50);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'labelMeLayer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'labelMeLayer1PhaseAll';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode1';
params.PhaseSimulation.weightAll = 200;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;

gridjobs = Gridjob(params);
start(gridjobs);
