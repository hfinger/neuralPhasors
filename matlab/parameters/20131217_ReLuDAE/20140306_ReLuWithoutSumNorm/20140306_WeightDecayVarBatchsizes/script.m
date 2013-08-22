clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ReLuDAE100';
params.Gridjob.initRandStreamWithSeed = 12345;

params.Autoencoder.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'ReLuDAE100';
params.Autoencoder.continue = false;
params.Autoencoder.debugOutput = true;

params.Autoencoder.continueBatch = false;
params.Autoencoder.continueBatchInWeightsFolder = [];
params.Autoencoder.continueBatchInBackConnFilenames = [];
params.Autoencoder.continueBatchInForwConnFilenames = [];

params.Autoencoder.inSamplesDims = {[60 60 3 200],[60 60 3 600],[60 60 3 2000]}; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 8;
params.Autoencoder.patchDimBackward = 8;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = 2;
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = false;
params.Autoencoder.useNormMean = false;
params.Autoencoder.useRectifiedLinear = true;
params.Autoencoder.saveinterval = 10;
params.Autoencoder.saveintervalpng = 10;
params.Autoencoder.savepng = true;
params.Autoencoder.sparsityParam = 0.01;
params.Autoencoder.alpha = 1;        % weight of Autoencoder
params.Autoencoder.lambda = 1;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
params.Autoencoder.beta = 0;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.gammaPvalue = 0.5;
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0.75;

params.Autoencoder.numberOfPatchReloads = 100;
params.Autoencoder.numberOfImagesPerPatchReload = 100;

params.Autoencoder.batchsize = []; %[] means full batch
params.Autoencoder.fixedBatches = false;
params.Autoencoder.validationSetsize = 0;
params.Autoencoder.validationInterval = Inf;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.validationSetImageIds = [];

params.Autoencoder.resetRandstreamEachReload = false;
params.Autoencoder.resetRandstreamEachIter = false;
params.Autoencoder.resetRandstreamEachEval = true;

params.Autoencoder.useMinFuncGrad = false;

params.Autoencoder.topoNghFcn = [];
params.Autoencoder.topoNgh = [1];
params.Autoencoder.topoGridDims = [100];
params.Autoencoder.topoPeriodicBoundary = [false];
params.Autoencoder.topoEpsilon = 1e-2;

params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = {5, 20, 50};
params.minFunc.display = 'on';
params.minFunc.optTol = 1e-6;
params.minFunc.progTol = 1e-14;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


