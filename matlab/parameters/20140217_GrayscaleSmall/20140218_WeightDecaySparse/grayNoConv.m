clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'simAnnealing';
params.Gridjob.initRandStreamWithSeed = 12345;

params.RBM.inActFolder = '../labelMeWhite';
params.RBM.inActFilenames = 'act.*.mat';
params.RBM.outWeightsFolder = 'RBMWeights'; %relative to the workpath
params.RBM.continue = false;

params.RBM.inSamplesDims = [12 12 1]; % [x,y,#features]
params.RBM.inSamplesBorderBuffer = 0;

params.RBM.patchDim = 1; % after subsampling of input
params.RBM.hiddenSize = 200;
params.RBM.inputSubsampling = 12;

params.RBM.initWMax = 0.01;
params.RBM.initWGaussian = true;
params.RBM.initWScaleDistFcn = [];

params.RBM.saveinterval = 10;
params.RBM.saveintervalmat = 10;
params.RBM.saveintervalpng = 10;

params.RBM.topo = 0;
params.RBM.topoNghFcn = [];
params.RBM.topoNgh = [3 3];
params.RBM.topoGridDims = [20 20];
params.RBM.topoPeriodicBoundary = [true true];
params.RBM.topoEpsilon = 1e-2;

params.RBM.numberOfPatchReloads = 10000;
params.RBM.numberOfImagesPerPatchReload = 100;
params.RBM.numberOfPatchesPerReload = 1000;
params.RBM.numberOfIterationsPerReload = 10;

params.RBM.validationSetsize = 0;
params.RBM.validationInterval = Inf;

params.RBM.lRate = 0.0003;			% LEARNING RATE
params.RBM.nGibbs = 1;				% # OF GIBBS SAMPLES (CONTRASTIVE DIVERGENCE)

params.RBM.sparsity = 0.01;		% TARGET HIDDEN UNIT SPARSITY
params.RBM.sparseGain = {0.003, 0.01, 0.03};			% GAIN ON THE LEARNING RATE FOR SPARSITY CONSTRAINTS
params.RBM.momentum = 0.9;			% (DEFAULT) GRADIENT MOMENTUM FOR WEIGHT UPDATES
params.RBM.wPenalty = {0.003, 0.01, 0.03};			% L2 WEIGHT PENALTY
params.RBM.beginAnneal = 400;		% BEGIN SIMULUATED ANNEALING AFTER THIS # OF EPOCHS
params.RBM.beginWeightDecay = 1;	% BEGIN WEIGHT DECAY AFTER THIS # OF EPOCHS


paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


