clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'lMLayer1noRectify';
params.ApplyWeights.inActFolder = 'labelMeWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.outActFolder = 'labelMeActLayer1noRectify';
params.ApplyWeights.weightsFile = 'cifarAutoencoderWeights/forwConn.mat';
% params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), min(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{1} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
