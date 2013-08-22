clear paramsAll;

clear params;

% params.Gridjob.runLocal = true;
% params.Gridjob.requiremf = 12000;
% params.Gridjob.jobname = 'rfs1WithoutActFcn';
% params.ApplyWeights.plotPdf = true;
% params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
% params.ApplyWeights.inActFilenames = 'act.*.mat';
% params.ApplyWeights.fileid = 1;
% params.ApplyWeights.outActFolder = 'rfs3WithoutActFcn';
% params.ApplyWeights.weightsFile = 'ReLuDAE100/1/forwConn.mat';
% params.ApplyWeights.convType = 'same';
% params.ApplyWeights.actFcn = @(x) x; 
% params.ApplyWeights.actFcn2 = @(x) x;
% paramsAll{1} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfs3WithActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../rfsWithActFcn/2';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfs3WithActFcn';
params.ApplyWeights.weightsFile = 'ReLuDAE200layer2/3/forwConn.mat';
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;

clear params;

gridjobs = Gridjob(paramsAll);
start(gridjobs);

