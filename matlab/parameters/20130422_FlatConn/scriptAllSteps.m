clear W;

rmax = 10;

dx = (-rmax:rmax)';
dy = (-rmax:rmax);

dx = repmat(dx,[1, size(dy,2)]);
dy = repmat(dy,[size(dx,1) 1]);

dx = dx(:);
dy = dy(:);

r = sqrt(dx.^2 + dy.^2);

W.dx = dx(r<=rmax);
W.dy = dy(r<=rmax);

W.f0 = ones(size(W.dx));
W.f1 = ones(size(W.dx));
W.w = ones(size(W.dx));

for i=1:10
  allW{i} = W;
  p = randperm(length(W.dx));
  useConnIds = p(1:round(i*length(W.dx)/10));
  allW{i}.dx = allW{i}.dx(useConnIds);
  allW{i}.dy = allW{i}.dy(useConnIds);
  allW{i}.f0 = allW{i}.f0(useConnIds);
  allW{i}.f1 = allW{i}.f1(useConnIds);
  allW{i}.w = allW{i}.w(useConnIds);
end

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1Phase';
params.PhaseSimulation.inActFolder = ones(1000,1000);
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = allW;
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
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
params.PhaseSimulation.maxdphase = 0.5;

gridjobs = Gridjob(params);
start(gridjobs);
