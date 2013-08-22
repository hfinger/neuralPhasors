clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'lMLayer1PhaseAllRK4';
params.Gridjob.wc_host = '';
params.PhaseSimulation.inActFolder = 'labelMeActLayer1';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(3:50);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'labelMeLayer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'labelMeLayer1PhaseAllRK4';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 1;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;

gridjobs = Gridjob(params);
start(gridjobs);
