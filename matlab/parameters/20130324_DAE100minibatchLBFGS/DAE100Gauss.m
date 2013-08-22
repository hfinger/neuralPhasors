clear params;

params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'DAE100Gauss';

params.Autoencoder.inActFolder = '../20130319_AE100LabelMe/labelMeWhite/05june05_static_street_boston/p1010736.jpg';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'DAE100Gauss';
params.Autoencoder.continue = false;

params.Autoencoder.inSamplesDims = [32 32 3 10000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;

params.Autoencoder.patchDimForward = 6;
params.Autoencoder.patchDimBackward = 4;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = 2;

params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];

params.Autoencoder.useSoftmax = true;
params.Autoencoder.useRectifiedLinear = false;%{false,true};

params.Autoencoder.saveinterval = 1;
params.Autoencoder.savepng = true;
params.Autoencoder.sparsityParam = 0.035;

params.Autoencoder.alpha = 1;        % weight of Autoencoder
params.Autoencoder.lambda = 0;    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0 ;
params.Autoencoder.beta = 0;         % weight of sparsity penalty term
params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = 0.5;

params.Autoencoder.batchsize = 500; %[] means full batch
params.Autoencoder.fixedBatches = false;

params.Autoencoder.validationSetsize = 500;
params.Autoencoder.validationInterval = 10;
params.Autoencoder.validationSetIds = [];

params.Autoencoder.useMinFuncGrad = false;

params.minFunc.Method = 'lbfgs';
params.minFunc.maxIter = 400;
params.minFunc.display = 'on';
params.minFunc.Corr = 5;
params.minFunc.Damped = 1;
params.minFunc.optTol = 0;
params.minFunc.progTol = -100;
params.minFunc.LS_interp = {0,1,2};
params.minFunc.LS_type = {0,1};

autoencoderJob = Autoencoder(params);
start(autoencoderJob);

