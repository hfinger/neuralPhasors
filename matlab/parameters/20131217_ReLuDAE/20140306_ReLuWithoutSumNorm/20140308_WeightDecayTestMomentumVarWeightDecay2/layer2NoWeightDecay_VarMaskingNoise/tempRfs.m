clear paramsAll;

clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfsWithoutActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../rfsWithActFcn/1';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfsWithoutActFcn';
params.ApplyWeights.weightsFile = {'temp_ReLuDAE200layer2/1/forwConnIter2000.mat','temp_ReLuDAE200layer2/2/forwConnIter2000.mat','temp_ReLuDAE200layer2/3/forwConnIter2000.mat','temp_ReLuDAE200layer2/4/forwConnIter2000.mat'};
params.ApplyWeights.convType = 'same';
params.ApplyWeights.actFcn = @(x) x; 
params.ApplyWeights.actFcn2 = @(x) x;
paramsAll{1} = params;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'rfsWithActFcn';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = '../rfsWithActFcn/1';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'rfsWithActFcn';
params.ApplyWeights.weightsFile = {'temp_ReLuDAE200layer2/1/forwConnIter2000.mat','temp_ReLuDAE200layer2/2/forwConnIter2000.mat','temp_ReLuDAE200layer2/3/forwConnIter2000.mat','temp_ReLuDAE200layer2/4/forwConnIter2000.mat'};
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;

gridjobs = Gridjob(paramsAll);
start(gridjobs);

