clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'AE100subsVarManySamples';
params.Gridjob.wc_host = [];

params.Autoencoder.inActFolder = 'labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'AE100subsVarManySamplesWeights';
params.Autoencoder.continue = false;

params.Autoencoder.inSamplesDims = [60 60 3 5000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;

params.Autoencoder.patchDimForward = 4;
params.Autoencoder.patchDimBackward = 4;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = 3;

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
params.Autoencoder.lambda = 0.003;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0 ;
params.Autoencoder.beta = 0.3;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0.2;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0;

params.Autoencoder.batchsize = 500; %[] means full batch
params.Autoencoder.fixedBatches = false;

params.Autoencoder.validationSetsize = 2500;
params.Autoencoder.validationInterval = 10;
params.Autoencoder.validationSetIds = [];

params.Autoencoder.resetRandstreamEachIter = true;
params.Autoencoder.useMinFuncGrad = false;

params.minFuncGrad.Method = 'rmsprop';
params.minFuncGrad.maxIter = 10000;
params.minFuncGrad.learnrate = 1e-3;
params.minFuncGrad.EMAconst = 0.1;
params.minFuncGrad.momentum = 0.9; % = 1 - 2/(lambda*2.8854+1) for halflife lambda
params.minFuncGrad.display = 'on';

params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 400;
params.minFunc.display = 'on';

autoencoderJob = Autoencoder(params);
start(autoencoderJob);

