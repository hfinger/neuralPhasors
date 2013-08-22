clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'layer2AE100weights';
params.Autoencoder.inActFolder = '../../20130726_Paper/Autoencoder/layer1ActNotRectified';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'layer2AE100weights';
params.Autoencoder.continue = false;
params.Autoencoder.inSamplesDims = [30 30 100 1000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 5;
params.Autoencoder.patchDimBackward = 5;
params.Autoencoder.hiddenSize = 200;
params.Autoencoder.inputSubsampling = 2;
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = false;
params.Autoencoder.useRectifiedLinear = false;
params.Autoencoder.reconstrSigmoid = true;
params.Autoencoder.useCrossEntropyError = true;
params.Autoencoder.saveinterval = 10;
params.Autoencoder.saveintervalpng = 10;
params.Autoencoder.savepng = false;
params.Autoencoder.sparsityParam = 0.035;
params.Autoencoder.alpha = 0.01/3;        % weight of Autoencoder
params.Autoencoder.lambda = {0.001, 0.003};    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
params.Autoencoder.beta = {0.1, 0.3};         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0;

params.Autoencoder.numberOfPatchReloads = 50;
params.Autoencoder.numberOfImagesPerPatchReload = 200;

params.Autoencoder.batchsize = 500; %[] means full batch
params.Autoencoder.fixedBatches = false;
params.Autoencoder.validationSetsize = 500;
params.Autoencoder.validationInterval = 10;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.resetRandstreamEachIter = true;
params.Autoencoder.useMinFuncGrad = false;
params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 20;
params.minFunc.display = 'on';
params.minFunc.optTol = 1e-6;
params.minFunc.progTol = 1e-14;
paramsAll{1} = params;


clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActNotRectified';
params.ApplyWeights.inActFolder = '../../20130726_Paper/Autoencoder/layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer2ActNotRectified';
params.ApplyWeights.weightsFile = {'temp_layer2AE100weights/1/forwConnIter270.mat','temp_layer2AE100weights/2/forwConnIter270.mat','temp_layer2AE100weights/3/forwConnIter270.mat','temp_layer2AE100weights/4/forwConnIter270.mat'};
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2ActRectified';
params.ApplyWeights.inActFolder = '../../20130726_Paper/Autoencoder/layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:185;
params.ApplyWeights.outActFolder = 'layer2ActRectified';
params.ApplyWeights.weightsFile = {'temp_layer2AE100weights/1/forwConnIter270.mat','temp_layer2AE100weights/2/forwConnIter270.mat','temp_layer2AE100weights/3/forwConnIter270.mat','temp_layer2AE100weights/4/forwConnIter270.mat'};
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
gridjobs = Gridjob(paramsAll(2:3));
start(gridjobs);


