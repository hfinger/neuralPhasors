clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeInput';
params.LoadLabelMe.catName = {'05june05_static_street_boston','05june05_static_street_porter'};
params.LoadLabelMe.fileid = [];
params.LoadLabelMe.outActFolder = 'labelMeInput';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeWhiteningWeights';
params.Gridjob.requiremf = 12000;
params.Whitening.inActFolder = 'labelMeInput';
params.Whitening.outWeightsFolder = 'labelMeWhiteningWeights';
params.Whitening.inNumChannels = 3;
params.Whitening.maxCorrLengthDim1 = 64;
params.Whitening.maxCorrLengthDim2 = 64;
params.Whitening.convKernelDim1 = 51;
params.Whitening.convKernelDim2 = 51;
params.Whitening.reduceToConv = true;
params.Whitening.epsilon = 0.1;
params.Whitening.numPatches = 10000;
params.Whitening.matchFilenames = 'act.*.mat';
params.Whitening.borderBuffer = 0;
paramsAll{2} = params;

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
paramsAll{3} = params;

clear params;
gridjobs = Gridjob(paramsAll(2:3));
start(gridjobs);


