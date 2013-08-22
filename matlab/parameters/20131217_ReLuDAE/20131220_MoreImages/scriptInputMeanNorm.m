clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeMeanNorm';
params.ApplyWeights.inActFolder = 'labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = [];
params.ApplyWeights.outActFolder = 'labelMeMeanNorm';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn = @(x) bsxfun(@rdivide,x,sum(x,3)+1e-6);
paramsAll{1} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
