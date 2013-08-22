clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'AE';

params.Autoencoder.inActFolder = 'cifarWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'cifarAutoencoderWeights';
params.Autoencoder.continue = false;

params.Autoencoder.inSamplesDims = [32 32 3 500]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;

params.Autoencoder.patchDimForward = 12;
params.Autoencoder.patchDimBackward = 12;
params.Autoencoder.hiddenSize = 400;
params.Autoencoder.inputSubsampling = 1;

params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];

params.Autoencoder.useSoftmax = false;
params.Autoencoder.useRectifiedLinear = false;

params.Autoencoder.saveinterval = 10;
params.Autoencoder.savepng = true;
params.Autoencoder.sparsityParam = 0.035;

params.Autoencoder.alpha = 0.01/3;        % weight of Autoencoder
params.Autoencoder.lambda = 0.001;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0 ;
params.Autoencoder.beta = 0.3;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0;

params.Autoencoder.topoNghFcn = [];
params.Autoencoder.topoNgh = [3 3];
params.Autoencoder.topoGridDims = [20 20];
params.Autoencoder.topoPeriodicBoundary = [true true];
params.Autoencoder.topoEpsilon = 1e-2;

params.Autoencoder.batchsize = []; %[] means full batch
params.Autoencoder.fixedBatches = false;

params.Autoencoder.validationSetsize = [];
params.Autoencoder.validationInterval = Inf;
params.Autoencoder.validationSetIds = [];

params.Autoencoder.useMinFuncGrad = false;

params.minFuncGrad.Method = 'rmsprop';
params.minFuncGrad.maxIter = 10000;
params.minFuncGrad.learnrate = 1e-5;
params.minFuncGrad.EMAconst = 0.1;
params.minFuncGrad.momentum = 0.9; % = 1 - 2/(lambda*2.8854+1) for halflife lambda
params.minFuncGrad.display = 'on';

params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 400;
params.minFunc.display = 'on';

autoencoderJob = Autoencoder(params);
start(autoencoderJob);

