clear W;

radmax = [2 3 6 10];
for i=1:4
  rmax = radmax(i);

  dx = (-rmax:rmax)';
  dy = (-rmax:rmax);

  dx = repmat(dx,[1, size(dy,2)]);
  dy = repmat(dy,[size(dx,1) 1]);

  dx = dx(:);
  dy = dy(:);

  r = sqrt(dx.^2 + dy.^2);

  allW{i}.dx = dx(r<=rmax);
  allW{i}.dy = dy(r<=rmax);

  allW{i}.f0 = ones(size(allW{i}.dx));
  allW{i}.f1 = ones(size(allW{i}.dx));
  allW{i}.w = ones(size(allW{i}.dx));

end

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1PhaseVarMaxRad5';
params.PhaseSimulation.inActFolder = ones(200,200);
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = allW;
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1PhaseVarMaxRad5';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 200;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 0.005;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;

gridjobs = Gridjob(params);
start(gridjobs);
