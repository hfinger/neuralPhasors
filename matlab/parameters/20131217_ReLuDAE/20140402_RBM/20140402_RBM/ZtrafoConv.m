clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'RBM';
params.Gridjob.initRandStreamWithSeed = 12345;

params.RBM.inActFolder = '../../20131220_MoreImages/labelMeWhite';
params.RBM.inActFilenames = 'act.*.mat';
params.RBM.outWeightsFolder = 'RBM'; %relative to the workpath
params.RBM.continue = false;

params.RBM.inSamplesDims = [40 40 3]; % [x,y,#features]
params.RBM.inSamplesBorderBuffer = 0;

params.RBM.patchDim = 4; % after subsampling of input
params.RBM.hiddenSize = 100;
params.RBM.inputSubsampling = 4;

params.RBM.initWMax = 0.01;
params.RBM.initWGaussian = true;
params.RBM.initWScaleDistFcn = [];

params.RBM.saveinterval = 50;
params.RBM.saveintervalmat = 50;
params.RBM.saveintervalpng = 50;

params.RBM.topo = 0;
params.RBM.topoNghFcn = [];
params.RBM.topoNgh = [3 3];
params.RBM.topoGridDims = [20 20];
params.RBM.topoPeriodicBoundary = [true true];
params.RBM.topoEpsilon = 1e-2;

params.RBM.numberOfPatchReloads = 10000;
params.RBM.numberOfImagesPerPatchReload = 100;
params.RBM.numberOfPatchesPerReload = 1000;
params.RBM.numberOfIterationsPerReload = 5;

params.RBM.validationSetsize = 0;
params.RBM.validationInterval = Inf;

params.RBM.lRate = 0.001;			% LEARNING RATE
params.RBM.nGibbs = 1;				% # OF GIBBS SAMPLES (CONTRASTIVE DIVERGENCE)

params.RBM.sparsity = 0.05;		% TARGET HIDDEN UNIT SPARSITY
params.RBM.sparseGain = {0.03, 0.1, 0.3};			% GAIN ON THE LEARNING RATE FOR SPARSITY CONSTRAINTS
params.RBM.momentum = 0.9;			% (DEFAULT) GRADIENT MOMENTUM FOR WEIGHT UPDATES
params.RBM.wPenalty = {0.01, 0.03, 0.1};%.05;			% L2 WEIGHT PENALTY
params.RBM.beginAnneal = 400;		% BEGIN SIMULUATED ANNEALING AFTER THIS # OF EPOCHS
params.RBM.beginWeightDecay = 1;	% BEGIN WEIGHT DECAY AFTER THIS # OF EPOCHS


paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


