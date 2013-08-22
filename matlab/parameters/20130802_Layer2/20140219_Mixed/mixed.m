clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'weightsL1toL2tmp';
params.CombineWeights.weightsFileIn1WithoutActFcn = '../NoPhaseFromLayer1/layer1ActRectifiedWhiteningWeights/weights.mat'; %relative to the workpath
params.CombineWeights.weightsFileIn2WithActFCN = '../NoPhaseFromLayer1/layer2AE100weights/2/forwConn.mat'; %relative to the workpath
params.CombineWeights.weightsFileOutCombined = 'weightsL1toL2tmp'; %relative to the workpath
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActNotRectified';
params.ApplyWeights.inActFolder = '../NoPhaseFromLayer1/layer1ActRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer2ActNotRectified';
params.ApplyWeights.weightsFile = 'weightsL1toL2/forwConn.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


