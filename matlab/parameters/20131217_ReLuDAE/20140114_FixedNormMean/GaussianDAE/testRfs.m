clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActWithoutActFcn';
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'layer1ActWithoutActFcn';
params.ApplyWeights.weightsFile = '../ReLuDAE100weightsVarNoise/4/forwConn.mat';
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn = @(x) x; 
params.ApplyWeights.actFcn2 = @(x) x;
paramsAll{1} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);