clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'testlayer1ActWithoutActFcn3';
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:10;
params.ApplyWeights.outActFolder = 'testlayer1ActWithoutActFcn3';
params.ApplyWeights.weightsFile = '../ReLuDAE100weightsSubs2SparseAct/3/forwConn.mat';
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn = @(x) x; 
params.ApplyWeights.actFcn2 = @(x) x;
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);

clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'testlayer1Act3';
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:10;
params.ApplyWeights.outActFolder = 'testlayer1Act3';
params.ApplyWeights.weightsFile = '../ReLuDAE100weightsSubs2SparseAct/3/forwConn.mat';
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);