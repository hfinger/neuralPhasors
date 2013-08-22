clear params;

params{1}.Gridjob.runLocal = true;
params{1}.Gridjob.jobname = 'lMload';
params{1}.LoadLabelMe.catName = '05june05_static_street_boston';
params{1}.LoadLabelMe.fileid = 1:100;
params{1}.LoadLabelMe.outActFolder = 'labelMeInput';
params{1}.LoadLabelMe.resizeToX = 200;
params{1}.LoadLabelMe.resizeToY = 150;

params{2}.Gridjob.runLocal = true;
params{2}.Gridjob.jobname = 'lMWhiteMat';
params{2}.Gridjob.requiremf = 12000;
params{2}.Whitening.inActFolder = 'labelMeInput';
params{2}.Whitening.outWeightsFolder = 'labelMeWhiteningWeights';
params{2}.Whitening.inNumChannels = 3;
params{2}.Whitening.maxCorrLengthDim1 = 64;
params{2}.Whitening.maxCorrLengthDim2 = 64;
params{2}.Whitening.convKernelDim1 = 51;
params{2}.Whitening.convKernelDim2 = 51;
params{2}.Whitening.reduceToConv = true;
params{2}.Whitening.epsilon = 0.1;
params{2}.Whitening.numPatches = 10000;
params{2}.Whitening.matchFilenames = 'act.*.mat';
params{2}.Whitening.borderBuffer = 0;

params{3}.Gridjob.runLocal = true;
params{3}.Gridjob.requiremf = 13000;
params{3}.Gridjob.jobname = 'lMWhitening';
params{3}.ApplyWeights.inActFolder = 'labelMeInput';
params{3}.ApplyWeights.inActFilenames = 'act.*.mat';
params{3}.ApplyWeights.fileid = num2cell(1:100);
params{3}.ApplyWeights.outActFolder = 'labelMeWhite';
params{3}.ApplyWeights.weightsFile = 'labelMeWhiteningWeights/weights.mat';

gridjobs = Gridjob(params);
start(gridjobs);


