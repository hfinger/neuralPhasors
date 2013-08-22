clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.jobname = 'labelMeWhiteningWeights';
params.Gridjob.requiremf = 12000;
params.Whitening.inActFolder = 'labelMeZtransformed';
params.Whitening.outWeightsFolder = 'labelMeWhiteningWeights';
params.Whitening.inNumChannels = 3;
params.Whitening.maxCorrLengthDim1 = 32;
params.Whitening.maxCorrLengthDim2 = 32;
params.Whitening.convKernelDim1 = 21;
params.Whitening.convKernelDim2 = 21;
params.Whitening.reduceToConv = true;
params.Whitening.epsilon = 0.1;
params.Whitening.numPatches = 10000;
params.Whitening.matchFilenames = 'act.*.mat';
params.Whitening.borderBuffer = 0;
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeWhite';
params.ApplyWeights.inActFolder = 'labelMeZtransformed';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = [];
params.ApplyWeights.outActFolder = 'labelMeWhite';
params.ApplyWeights.weightsFile = 'labelMeWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.jobname = 'labelMeWhiteZtrafoWeights';
params.Gridjob.requiremf = 12000;
params.PreprocessZtrafo.inActFolder = 'labelMeWhite';
params.PreprocessZtrafo.outWeightsFolder = 'labelMeWhiteZtrafoWeights';
params.PreprocessZtrafo.inNumChannels = 3;
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeWhiteZtransformed';
params.ApplyWeights.inActFolder = 'labelMeWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = [];
params.ApplyWeights.outActFolder = 'labelMeWhiteZtransformed';
params.ApplyWeights.weightsFile = 'labelMeWhiteZtrafoWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeConcat';
params.ActConcat.inActFolder1 = 'labelMeZtransformed';
params.ActConcat.inActFilenames1 = 'act.*.mat';
params.ActConcat.inActFolder2 = 'labelMeWhiteZtransformed';
params.ActConcat.inActFilenames2 = 'act.*.mat';
params.ActConcat.fileid = [];
params.ActConcat.outActFolder = 'labelMeConcat';
paramsAll{5} = params;

clear params;
gridjobs = Gridjob(paramsAll(3:5));
start(gridjobs);


