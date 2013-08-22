clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'PreprocessZtrafo';
params.Gridjob.requiremf = 12000;
params.PreprocessZtrafo.inActFolder = '../20140129_InputSumNormalized/labelMeInput';
params.PreprocessZtrafo.outWeightsFolder = 'zTrafoWeights';
params.PreprocessZtrafo.inNumChannels = 3;
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeZtransformed';
params.ApplyWeights.inActFolder = '../20140129_InputSumNormalized/labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = [];
params.ApplyWeights.outActFolder = 'labelMeZtransformed';
params.ApplyWeights.weightsFile = 'zTrafoWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


