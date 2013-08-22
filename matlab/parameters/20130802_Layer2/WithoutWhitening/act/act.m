clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActNotRectified';
params.ApplyWeights.inActFolder = '../../NoPhaseFromLayer1/layer1ActRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'layer2ActNotRectified';
params.ApplyWeights.weightsFile = '../temp_layer2AE100weights/forwConnIter300.mat';
params.ApplyWeights.convType = 'same';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.plotColormap = 'hot';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActRectified';
params.ApplyWeights.inActFolder = '../../NoPhaseFromLayer1/layer1ActRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1;
params.ApplyWeights.outActFolder = 'layer2ActRectified';
params.ApplyWeights.weightsFile = '../temp_layer2AE100weights/forwConnIter300.mat';
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
params.ApplyWeights.plotPdf = true;
params.ApplyWeights.plotColormap = 'hot';
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
