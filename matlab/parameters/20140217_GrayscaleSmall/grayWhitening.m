clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeWhiteningWeights';
params.Gridjob.requiremf = 12000;
params.Whitening.inActFolder = 'labelMeInput';
params.Whitening.outWeightsFolder = 'labelMeWhiteningWeights';
params.Whitening.inNumChannels = 1;
params.Whitening.maxCorrLengthDim1 = 32;
params.Whitening.maxCorrLengthDim2 = 32;
params.Whitening.convKernelDim1 = 12;
params.Whitening.convKernelDim2 = 12;
params.Whitening.reduceToConv = true;
params.Whitening.epsilon = 0.0001;
params.Whitening.numPatches = 10000;
params.Whitening.matchFilenames = 'act.*.mat';
params.Whitening.borderBuffer = 0;
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeWhite';
params.ApplyWeights.inActFolder = 'labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = [];
params.ApplyWeights.outActFolder = 'labelMeWhite';
params.ApplyWeights.weightsFile = 'labelMeWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;


clear params;
gridjobs = Gridjob(paramsAll(2));
start(gridjobs);


