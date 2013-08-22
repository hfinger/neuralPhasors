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
conn.W = cat(6,W,-W); % add opposite colors
conn.inputSubtract = zeros(1,1,3);
conn.inputSubtract(1,1,:) = [0.4045, 0.4165, 0.4027];

%%
clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'lMLayer1';
params.ApplyWeights.inActFolder = '../20130318_AE400/labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'labelMeActLayer1';
params.ApplyWeights.weightsFile = conn;%'forwConn.mat';
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.jobname = 'lMLayer1Cov';
params.Gridjob.requiremf = 13000;
params.FeatureCovariance.inActFolder = 'labelMeActLayer1';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = [];
params.FeatureCovariance.outCovFolder = 'labelMeLayer1Cov';
params.FeatureCovariance.maxCovLengthDim1 = 32;
params.FeatureCovariance.maxCovLengthDim2 = 32;
params.FeatureCovariance.numSamplesPerImage = 20;
params.FeatureCovariance.borderBuffer = 0;
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.jobname = 'lMLayer1Conn';
params.Gridjob.requiremf = 13000;
params.ProbabilisticConn.inCorrFile = 'labelMeLayer1Cov/patchCorr.mat';
params.ProbabilisticConn.outWeightsFolder = 'labelMeLayer1Conn';
params.ProbabilisticConn.numExc = 1000;
params.ProbabilisticConn.numInh = 200;
params.ProbabilisticConn.exclSelfConn = true;
params.ProbabilisticConn.useDiscretesample = true;
params.ProbabilisticConn.exclDoubleConns = true;
params.ProbabilisticConn.probRadiusFcn = [];
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'lMLayer1Phase';
params.PhaseSimulation.inActFolder = 'labelMeActLayer1';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = 1;
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = 'labelMeLayer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'labelMeLayer1Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode1';
params.PhaseSimulation.weightAll = {0.1,1,10,100};
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 1;
params.PhaseSimulation.saveintervalMeanPhase = 1;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;
paramsAll{4} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
