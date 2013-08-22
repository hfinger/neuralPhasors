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
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'lMLayer1noRectify';
params.ApplyWeights.inActFolder = '../20130318_AE400/labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'labelMeActLayer1noRectify';
params.ApplyWeights.weightsFile = conn;%'forwConn.mat';
params.ApplyWeights.actFcn = @(x) 1 ./ (1 + exp(-x));
% params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), min(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
