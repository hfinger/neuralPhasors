clear paramsAll;

clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfs3WithoutActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfs3WithoutActFcn';
params.ApplyWeights.weightsFile = 'temp_ReLuDAE100/3/forwConnIter650.mat';
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn = @(x) x; 
params.ApplyWeights.actFcn2 = @(x) x;
paramsAll{1} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfs3WithActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfs3WithActFcn';
params.ApplyWeights.weightsFile = 'temp_ReLuDAE100/3/forwConnIter650.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;

gridjobs = Gridjob(paramsAll);
start(gridjobs);

