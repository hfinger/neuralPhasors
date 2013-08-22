clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'lMLayer1PhaseGlobInh4';
params.PhaseSimulation.inActFolder = 'labelMeActLayer1';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'labelMeLayer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'labelMeLayer1PhaseGlobInh4';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 25;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode1';
params.PhaseSimulation.weightAll = 1; %without weightGlobInh
params.PhaseSimulation.weightInh = 0;
params.PhaseSimulation.weightExc = 100;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;
params.PhaseSimulation.weightInhRadialFcn = {@(r) max(10-r,0), @(r) 1./max(2,r), @(r) 1./max(2,sqrt(r)), @(r) exp(-r.^2/(2*30^2))};
params.PhaseSimulation.weightGlobInh = 50; % 50 corresponds to equal to exc

gridjobs = Gridjob(params);
start(gridjobs);
