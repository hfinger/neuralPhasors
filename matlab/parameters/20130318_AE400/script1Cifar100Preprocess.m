clear params;

params{1}.Gridjob.runLocal = true;
params{1}.Gridjob.jobname = 'cifarLoad';
params{1}.LoadCifar100.outActFolder = 'cifarInput';

params{2}.Gridjob.runLocal = true;
params{2}.Gridjob.jobname = 'cifarWhiteMat';
params{2}.Gridjob.requiremf = 8000;
params{2}.Whitening.inActFolder = 'cifarInput';
params{2}.Whitening.outWeightsFolder = 'cifarWhiteningWeights';
params{2}.Whitening.inNumChannels = 3;
params{2}.Whitening.maxCorrLengthDim1 = 32;
params{2}.Whitening.maxCorrLengthDim2 = 32;
params{2}.Whitening.convKernelDim1 = 21;
params{2}.Whitening.convKernelDim2 = 21;
params{2}.Whitening.reduceToConv = false;
params{2}.Whitening.epsilon = 0.1;
params{2}.Whitening.numPatches = 500;
params{2}.Whitening.matchFilenames = 'act.*.mat';
params{2}.Whitening.borderBuffer = 0;

params{3}.Gridjob.runLocal = true;
params{3}.Gridjob.requiremf = 4000;
params{3}.Gridjob.jobname = 'cifarWhitening';
params{3}.ApplyWeights.inActFolder = 'cifarInput';
params{3}.ApplyWeights.inActFilenames = 'act.*.mat';
params{3}.ApplyWeights.outActFolder = 'cifarWhite';
params{3}.ApplyWeights.weightsFile = 'cifarWhiteningWeights/weights.mat';

gridjobs = Gridjob(params);
start(gridjobs);
