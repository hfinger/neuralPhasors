clear paramsAll;

folderlist = { ...
'05june05_static_street_boston',...
'05june05_static_street_porter',...
'10feb04_static_cars_highland',...
'30may05_static_street_cambridge',...
'april21_static_outdoor_davis',...
'april21_static_outdoor_kendall',...
'barcelona_static_street',...
'boston_static_march',...
'dec_static_street',...
'madrid_static_street',...
'nov25_static_street',...
'nov6_static_outdoor',...
'nov7_static_outdoor',...
'oct6_static_outdoor',...
'paris_static_street',...
'static_barcelona_street_city_outdoor_2005',...
'static_harvard_outdoor_street',...
'static_outdoor_anchorage_alaska_usa'};
% Removed: 
% static_outdoor_city_laredo_spain
% static_dartmouth_hanover_june_2006
% static_nature_web_outdoor_animal
% static_newyork_city_urban
% static_outdoor_bay_area_submitted_alyosha_efros
% static_outdoor_bozeman_montana_usa
% static_outdoor_city_laredo_spain


clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeInput';
params.LoadLabelMe.catName = folderlist;
params.LoadLabelMe.fileid = [];
params.LoadLabelMe.outActFolder = 'labelMeInput';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeWhiteningWeights';
params.Gridjob.requiremf = 12000;
params.Whitening.inActFolder = 'labelMeInput';
params.Whitening.outWeightsFolder = 'labelMeWhiteningWeights';
params.Whitening.inNumChannels = 3;
params.Whitening.maxCorrLengthDim1 = 64;
params.Whitening.maxCorrLengthDim2 = 64;
params.Whitening.convKernelDim1 = 51;
params.Whitening.convKernelDim2 = 51;
params.Whitening.reduceToConv = true;
params.Whitening.epsilon = 0.1;
params.Whitening.numPatches = 10000;
params.Whitening.matchFilenames = 'act.*.mat';
params.Whitening.borderBuffer = 0;
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeWhite';
params.ApplyWeights.inActFolder = 'labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = [];
params.ApplyWeights.outActFolder = 'labelMeWhite';
params.ApplyWeights.weightsFile = 'labelMeWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;


clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'ReLuDAE100weights';
params.Gridjob.initRandStreamWithSeed = 12345;

params.Autoencoder.inActFolder = 'labelMeWhite';
params.Autoencoder.inActFilenames = 'act.*.mat';
params.Autoencoder.outWeightsFolder = 'ReLuDAE100weights';
params.Autoencoder.continue = false;
params.Autoencoder.debugOutput = true;

params.Autoencoder.continueBatch = false;
params.Autoencoder.continueBatchInWeightsFolder = [];
params.Autoencoder.continueBatchInBackConnFilenames = [];
params.Autoencoder.continueBatchInForwConnFilenames = [];

params.Autoencoder.inSamplesDims = [60 60 3 2000]; % [x,y,#features,#samples]
params.Autoencoder.inSamplesBorderBuffer = 0;
params.Autoencoder.patchDimForward = 4;
params.Autoencoder.patchDimBackward = 4;
params.Autoencoder.hiddenSize = 100;
params.Autoencoder.inputSubsampling = 4;
params.Autoencoder.forwInitMaxWeight = 0.01;
params.Autoencoder.forwInitScaleDistFcn = [];
params.Autoencoder.backInitMaxWeight = 0.01;
params.Autoencoder.backInitScaleDistFcn = [];
params.Autoencoder.useSoftmax = true;
params.Autoencoder.useRectifiedLinear = true;
params.Autoencoder.saveinterval = 10;
params.Autoencoder.saveintervalpng = 1;
params.Autoencoder.savepng = true;
params.Autoencoder.sparsityParam = 0;
params.Autoencoder.alpha = 1;        % weight of Autoencoder
params.Autoencoder.lambda = {0, 0.03, 0.1, 0.3};    % weight of L2-decay parameter
params.Autoencoder.lambdaL1 = 0;    % weight of L1-decay parameter
params.Autoencoder.scaleLambdaL1WithDistExponent = 0;
params.Autoencoder.beta = 0;         % weight of sparsity penalty term
params.Autoencoder.gamma = {0, 0.1};        % weight of topographic sparsity penalty term
params.Autoencoder.maskingNoiseFraction = 0.25;        % fraction of inputs to set to 0
params.Autoencoder.gaussianNoiseSigmaInput = {0,0.5};

params.Autoencoder.numberOfPatchReloads = 50;
params.Autoencoder.numberOfImagesPerPatchReload = 200;

params.Autoencoder.batchsize = []; %[] means full batch
params.Autoencoder.fixedBatches = false;
params.Autoencoder.validationSetsize = 1000;
params.Autoencoder.validationInterval = 10;
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
params.minFunc.maxIter = 20;
params.minFunc.display = 'on';
params.minFunc.optTol = 1e-6;
params.minFunc.progTol = 1e-14;
paramsAll{4} = params;

clear params;
gridjobs = Gridjob(paramsAll(4));
start(gridjobs);


