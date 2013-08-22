clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'layer1ActRectifiedWhiteningWeights';
params.Gridjob.requiremf = 12000;
params.Gridjob.wc_host = '';
params.Whitening.inActFolder = 'layer1ActRectified';
params.Whitening.outWeightsFolder = 'layer1ActRectifiedWhiteningWeights';
params.Whitening.inNumChannels = 100;
params.Whitening.maxCorrLengthDim1 = 10;
params.Whitening.maxCorrLengthDim2 = 10;
params.Whitening.convKernelDim1 = 7;
params.Whitening.convKernelDim2 = 7;
params.Whitening.reduceToConv = true;
params.Whitening.epsilon = 0.1;
params.Whitening.numPatches = 1000;
params.Whitening.matchFilenames = 'act.*.mat';
params.Whitening.borderBuffer = 0;
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'layer1ActRectifiedWhite';
params.Gridjob.wc_host = '';
params.ApplyWeights.inActFolder = 'layer1ActRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = num2cell(1:185);
params.ApplyWeights.outActFolder = 'layer1ActRectifiedWhite';
params.ApplyWeights.weightsFile = 'layer1ActRectifiedWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 9000;
params.Gridjob.jobname = 'layer2AE100weights';
params.Gridjob.wc_host = '';
params.Autoencoder.inActFolder = 'layer1ActRectified';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'layer2AE100weights';
params.Autoencoder.continue = false;
params.Autoencoder.inSamplesDims = [30 30 100 2000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 6;
params.Autoencoder.patchDimBackward = 6;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = [];
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = false;
params.Autoencoder.useRectifiedLinear = false;
params.Autoencoder.saveinterval = 10;
params.Autoencoder.saveintervalpng = Inf;
params.Autoencoder.savepng = false;
params.Autoencoder.sparsityParam = 0.035;
params.Autoencoder.alpha = 0.01/3;        % weight of Autoencoder
params.Autoencoder.lambda = 0.001;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
params.Autoencoder.beta = 0.3;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0;
params.Autoencoder.batchsize = []; %[] means full batch
params.Autoencoder.fixedBatches = false;
params.Autoencoder.validationSetsize = 1000;
params.Autoencoder.validationInterval = 10;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.resetRandstreamEachIter = true;
params.Autoencoder.useMinFuncGrad = false;
params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 500;
params.minFunc.display = 'on';
params.minFunc.optTol = 1e-6;
params.minFunc.progTol = 1e-14;
paramsAll{3} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


