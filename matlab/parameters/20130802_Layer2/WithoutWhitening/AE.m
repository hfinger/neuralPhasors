clear paramsAll;


clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 9000;
params.Gridjob.jobname = 'layer2AE100weights';
params.Gridjob.wc_host = 'elara';
params.Autoencoder.inActFolder = '../NoPhaseFromLayer1/layer1ActRectified';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'layer2AE100weights';
params.Autoencoder.continue = false;
params.Autoencoder.inSamplesDims = [30 30 100 4000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 6;
params.Autoencoder.patchDimBackward = 6;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = 1;
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = false;
params.Autoencoder.useRectifiedLinear = false;
params.Autoencoder.saveinterval = 5;
params.Autoencoder.saveintervalpng = 5;
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
params.Autoencoder.batchsize = 1000; %[] means full batch
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
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


