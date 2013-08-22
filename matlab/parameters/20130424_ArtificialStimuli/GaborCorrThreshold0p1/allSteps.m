%% create Gabor filter weights
rfSize = 12;
numOrient = 8;
order = 2;
sigma = 1;
aspectratio = 3;
integral = 0.01;
numChan = 3;

W = zeros(rfSize,rfSize,numChan,1,1,numChan*numOrient);
Wtmp = zeros(rfSize,rfSize,1,1,1,numOrient);
orientations=(0:numOrient-1)*pi/numOrient;
for k=1:length(orientations);
  Wtmp(:,:,1,1,1,k) = genSpatialKernel( rfSize, order, sigma, aspectratio, integral, orientations(k) );
end
for k=1:numChan
  W(:,:,k,1,1,(k-1)*numOrient+1:k*numOrient) = Wtmp;
end
conn.W = 5 * cat(6,W,-W); % add opposite colors
% so W in total consists of surround -5, center 10.01 and surround -5
conn.inputSubtract = zeros(1,1,3);
conn.inputSubtract(1,1,:) = [0.4255, 0.4409, 0.4324]; %subtract mean
conn.inputScaling = zeros(1,1,3);
conn.inputScaling(1,1,:) = [1/0.2465, 1/0.2432, 1/0.2563]; %standardize input
conn.inputSubsampling = 2;

%%
clear paramsAll;
clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'imageInput';
params.Gridjob.wc_host = '';
params.LoadImage.filepath = {...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3pt.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3ptBroken.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3ptColorGradientOnGray.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3ptColorGradient.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/cube.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/Kanizsa_triangle.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/Kanizsa_triangle_only_dots.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/Kanizsa_triangle_only_triangle.png'};
params.LoadImage.outActFolder = 'imageInput';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1ActNotRectified';
params.ApplyWeights.inActFolder = 'imageInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'layer1ActNotRectified';
params.ApplyWeights.weightsFile = conn;
params.ApplyWeights.actFcn = @(x) 1 ./ (1 + exp(-x));
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1ActRectified';
params.ApplyWeights.inActFolder = 'layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), min(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'layer1Phase';
params.Gridjob.wc_host = '';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:9);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = '../../20130409_GaborNoWhiteningWithSigmoidWithZtrafo5timesWithSubsampling/excludeSmallCorrs/layer1Conn0p1/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 50;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 3;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;
paramsAll{4} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


