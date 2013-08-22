act = 0.6*ones(200,200,3);
rad = 70;
tmp=[];
for x=-99.5:99.5
  for y=-99.5:99.5
    r = sqrt(x^2+y^2);
    angle = atan2(y,x);
    tmp(end+1)=angle;
    if mod(angle*20/pi,2)>1
      if r < 73 && r > 67
        act(100.5+x,100.5+y,:) = max(0,min(0.6,0.3*(abs(r-70)-0.5)));
      end
    end
  end
end
imagefolder = '/net/store/nbp/phasesim/workdir/20130726_Paper/ArtificialStimuli/imageInput';
mkdir(imagefolder)
save(fullfile(imagefolder,'act.mat'),'act')

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
% conn.inputSubtract = zeros(1,1,3);
% conn.inputSubtract(1,1,:) = [0.4255, 0.4409, 0.4324]; %subtract mean
% conn.inputScaling = zeros(1,1,3);
% conn.inputScaling(1,1,:) = [1/0.2465, 1/0.2432, 1/0.2563]; %standardize input
conn.inputSubsampling = 2;

%mkdir('/net/store/nbp/phasesim/workdir/20130726_Paper/Gabor/GaborWeights')
%save('/net/store/nbp/phasesim/workdir/20130726_Paper/Gabor/GaborWeights/forwConn.mat','-struct','conn')

%%
clear paramsAll;


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
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1ActRectified';
params.ApplyWeights.inActFolder = 'layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1PhaseAllManySamples';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = '../Gabor/layer1ConnManySamples/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1PhaseAllManySamples';
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
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
paramsAll{3} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
