clear W;

rmax = 400;

dx = (-rmax:rmax)';
dy = (-rmax:rmax);

dx = repmat(dx,[1, size(dy,2)]);
dy = repmat(dy,[size(dx,1) 1]);

dx = dx(:);
dy = dy(:);

r = sqrt(dx.^2 + dy.^2);
r(r==0) = 1;

for i=1:5
  
  sigma = i*2;
  rayleigh = r/sigma^2 .* exp( - r.^2 / (2*sigma^2) );
  
  rayleighNormByR = rayleigh./r;
  
  connProb = rayleighNormByR;%1/(r.^i);
  connProb = connProb(:)/sum(connProb(:));
  
  connIds = discretesample(connProb, 200);
  
  W.dx = dx(connIds);
  W.dy = dy(connIds);

  W.f0 = ones(size(W.dx));
  W.f1 = ones(size(W.dx));
  W.w = ones(size(W.dx));

  allW{i} = W;
end

%%
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'simVarConnPropVsRadExpNorm';
params.PhaseSimulation.inActFolder = ones(1000,1000);
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = allW;
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'simVarConnPropVsRadExNorm';
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
