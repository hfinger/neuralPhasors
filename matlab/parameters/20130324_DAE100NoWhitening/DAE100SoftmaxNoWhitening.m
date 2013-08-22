clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'DAE100SoftmaxNoWhitening';

params.Autoencoder.inActFolder = 'labelMeMeanSubtract';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'DAE100SoftmaxNoWhitening';
params.Autoencoder.continue = false;

params.Autoencoder.inSamplesDims = [64 64 3 10000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;

params.Autoencoder.patchDimForward = 6;
params.Autoencoder.patchDimBackward = 2;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = 2;

params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];

params.Autoencoder.useSoftmax = true;

params.Autoencoder.useRectifiedLinear = {false,true};
params.Autoencoder.useNoActFcn = false;

params.Autoencoder.saveinterval = 10;
params.Autoencoder.savepng = false;
params.Autoencoder.sparsityParam = 0.035;

params.Autoencoder.alpha = 1;        % weight of Autoencoder
params.Autoencoder.lambda = 0;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0 ;
params.Autoencoder.beta = 0;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = {0, 0.25};        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0.5;

params.Autoencoder.batchsize = 500; %[] means full batch
params.Autoencoder.fixedBatches = false;

params.Autoencoder.validationSetsize = 500;
params.Autoencoder.validationInterval = 10;
params.Autoencoder.validationSetIds = [];

params.Autoencoder.useMinFuncGrad = true;

params.minFuncGrad.Method = 'rmsprop';
params.minFuncGrad.maxIter = 10000;
params.minFuncGrad.learnrate = 1e-2;
params.minFuncGrad.EMAconst = 0.1;
params.minFuncGrad.momentum = 0.9; % = 1 - 2/(lambda*2.8854+1) for halflife lambda
params.minFuncGrad.display = 'on';

params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 400;
params.minFunc.display = 'on';

autoencoderJob = Autoencoder(params);
start(autoencoderJob);

