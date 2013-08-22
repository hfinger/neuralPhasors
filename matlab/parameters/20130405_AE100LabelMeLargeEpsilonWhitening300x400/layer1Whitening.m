clear params;

params{2}.Gridjob.runLocal = false;
params{2}.Gridjob.jobname = 'layer1White';
params{2}.Gridjob.requiremf = 25000;
params{2}.Gridjob.wc_host = 'isonoe';
params{2}.Whitening.inActFolder = 'layer1Act';
params{2}.Whitening.outWeightsFolder = 'layer1WhiteningWeights';
params{2}.Whitening.inNumChannels = 100;
params{2}.Whitening.maxCorrLengthDim1 = 15;
params{2}.Whitening.maxCorrLengthDim2 = 15;
params{2}.Whitening.convKernelDim1 = 12;
params{2}.Whitening.convKernelDim2 = 12;
params{2}.Whitening.reduceToConv = true;
params{2}.Whitening.epsilon = 0.1;
params{2}.Whitening.numPatches = 10000;
params{2}.Whitening.matchFilenames = 'act.*.mat';
params{2}.Whitening.borderBuffer = 0;

params{3}.Gridjob.runLocal = false;
params{3}.Gridjob.requiremf = 13000;
params{3}.Gridjob.jobname = 'layer1Whitening';
params{3}.ApplyWeights.inActFolder = 'layer1Act';
params{3}.ApplyWeights.inActFilenames = 'act.*.mat';
params{3}.ApplyWeights.fileid = 1:185;
params{3}.ApplyWeights.outActFolder = 'layer1White';
params{3}.ApplyWeights.weightsFile = 'layer1WhiteningWeights/weights.mat';
params{3}.ApplyWeights.convType = 'same';

% params(1) = [];
gridjobs = Gridjob(params(2:3));
start(gridjobs);


