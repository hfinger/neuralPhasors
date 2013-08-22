%% create Gabor filter weights
rfSize = 12;
numOrient = 8;
order = 2;
sigma = 1.5;
aspectratio = 2;
integral = 0.1;
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
% so W in total consists of surround -5, center 10.5 and surround -5
conn.inputSubtract = zeros(1,1,3);
conn.inputSubtract(1,1,:) = [0.5, 0.5, 0.5]; %subtract mean
% conn.inputScaling = zeros(1,1,3);
% conn.inputScaling(1,1,:) = [1/0.2465, 1/0.2432, 1/0.2563]; %standardize input
conn.inputSubsampling = 2;

%%
clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'imageInput';
params.Gridjob.wc_host = '';
params.LoadImage.filepath = {...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/circle2.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/circle3.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/circle4.png'};
params.LoadImage.outActFolder = 'imageInput';
paramsAll{1} = params;

% clear params;
% params.Gridjob.runLocal = true;
% params.Gridjob.requiremf = 13000;
% params.Gridjob.jobname = 'imagesWhite';
% params.ApplyWeights.inActFolder = 'imageInput';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.fileid = 1:2;
% params.ApplyWeights.outActFolder = 'imagesWhite';
% params.ApplyWeights.weightsFile = '../../20130409_AE100LabelMeLargeEpsilonWhiteningWithSubsampling/labelMeWhiteningWeights/weights.mat';
% params.ApplyWeights.convType = 'same';
% paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1ActNotRectified';
params.Gridjob.wc_host = '';
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
params.Gridjob.wc_host = '';
params.ApplyWeights.inActFolder = 'layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1PhaseAllManySamples3';
params.Gridjob.wc_host = '';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:3);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = '../layer1ConnManySamples/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1PhaseAllManySamples3';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 300;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 3;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
paramsAll{4} = params;

gridjobs = Gridjob(paramsAll);
start(gridjobs);
