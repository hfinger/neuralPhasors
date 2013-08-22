clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'testlayer1ActWithoutActFcn3';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = 'labelMeMeanNorm';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'testlayer1ActWithoutActFcn';
params.ApplyWeights.weightsFile = cellfun(@(x) ['temp_ReLuDAE100weightsSubs2Patch8VarNoiseVarGamma/' x '/forwConn.mat'],num2cell(1:3),'UniformOutput',false);
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
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.inActFolder = 'labelMeMeanNorm';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'testlayer1Act';
params.ApplyWeights.weightsFile = cellfun(@(x) ['temp_ReLuDAE100weightsSubs2Patch8VarNoiseVarGamma/' x '/forwConn.mat'],num2cell(1:3),'UniformOutput',false);
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);