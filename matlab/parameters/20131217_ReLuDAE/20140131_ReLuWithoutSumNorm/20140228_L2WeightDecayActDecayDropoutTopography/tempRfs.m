clear paramsAll;

clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfsWithoutActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfsWithoutActFcn';
params.ApplyWeights.weightsFile = cellfun(@(x) ['temp_ReLuDAE100/' num2str(x) '/forwConnIter340.mat'],num2cell([1 3]),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn = @(x) x; 
params.ApplyWeights.actFcn2 = @(x) x;
paramsAll{1} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfsWithActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../../20131220_MoreImages/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfsWithActFcn';
params.ApplyWeights.weightsFile = cellfun(@(x) ['temp_ReLuDAE100/' num2str(x) '/forwConnIter340.mat'],num2cell([1 3]),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;

gridjobs = Gridjob(paramsAll);
start(gridjobs);

