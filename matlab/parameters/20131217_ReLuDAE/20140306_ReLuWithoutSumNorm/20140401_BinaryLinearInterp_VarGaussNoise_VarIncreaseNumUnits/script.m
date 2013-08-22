clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
% params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'ReLuDAE100';
params.Gridjob.initRandStreamWithSeed = 12345;

params.Autoencoder.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'ReLuDAE100';
params.Autoencoder.continue = true;
params.Autoencoder.debugOutput = true;

params.Autoencoder.continueBatch = false;
params.Autoencoder.continueBatchInWeightsFolder = [];
params.Autoencoder.continueBatchInBackConnFilenames = [];
params.Autoencoder.continueBatchInForwConnFilenames = [];

params.Autoencoder.inSamplesDims = [60 60 3 2200]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 8;
params.Autoencoder.patchDimBackward = 8;
params.Autoencoder.hiddenSize = {100 ,ceil(100*(1:10).^0.5/10^0.5), ceil(100*(1:10).^1/10^1), ceil(100*(1:10).^2/10^2)};

plot(1:10,cell2mat(params.Autoencoder.hiddenSize(2:end)'))
legend('0.5','1','2')

params.Autoencoder.inputSubsampling = 2;
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = false;
params.Autoencoder.useNormMean = false;
params.Autoencoder.useRectifiedLinear = false;

params.Autoencoder.useBinaryLinearInterp = true;
params.Autoencoder.binaryLinearInterpScaling = 1;
params.Autoencoder.binaryLinearInterpEMAconst = 0.9;

params.Autoencoder.saveinterval = Inf;
params.Autoencoder.saveintervalForwConn = 500;
params.Autoencoder.saveintervalBackConn = Inf;
params.Autoencoder.saveintervalTrainInfo = Inf;
params.Autoencoder.saveintervalLog = 50;
params.Autoencoder.saveintervalpng = 100;
params.Autoencoder.reloadSaveinterval = 5;
params.Autoencoder.savepng = true;
params.Autoencoder.sparsityParam = 0.01;
params.Autoencoder.alpha = 1;        % weight of Autoencoder
params.Autoencoder.lambda = 0;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
params.Autoencoder.beta = 0;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.gammaPvalue = 0.5;
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = {0.3, 0.75, 1.5};

params.Autoencoder.numberOfPatchReloads = 20;
params.Autoencoder.numberOfImagesPerPatchReload = 300;

params.Autoencoder.batchsize = 10; %[] means full batch
params.Autoencoder.fixedBatches = false;
params.Autoencoder.validationSetsize = 200;
params.Autoencoder.validationInterval = Inf;
params.Autoencoder.validationIntervalLog = 50;
params.Autoencoder.validationSetIds = [];
params.Autoencoder.validationSetImageIds = 1:10:1001;

params.Autoencoder.resetRandstreamEachReload = false;
params.Autoencoder.resetRandstreamEachIter = false;
params.Autoencoder.resetRandstreamEachEval = false;

params.Autoencoder.useMinFuncGrad = true;

params.Autoencoder.topoNghFcn = [];
params.Autoencoder.topoNgh = [1];
params.Autoencoder.topoGridDims = [100];
params.Autoencoder.topoPeriodicBoundary = [false];
params.Autoencoder.topoEpsilon = 1e-2;

params.minFuncGrad.Method = 'momentum';
params.minFuncGrad.maxIter = 250;
params.minFuncGrad.learnrate = 3e-4;
params.minFuncGrad.momentum = 0.9;
params.minFuncGrad.display = 'on';
params.minFuncGrad.displayEvery = 10;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


