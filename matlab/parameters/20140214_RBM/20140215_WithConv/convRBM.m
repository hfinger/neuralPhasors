clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'cRBM100';
params.Gridjob.initRandStreamWithSeed = 12345;

params.RBM.inActFolder = '../../20131217_ReLuDAE/20140204_InputZtrafo/labelMeZtransformed';
params.RBM.inActFilenames = 'act.*.mat';
params.RBM.outWeightsFolder = 'cRBMWeights'; %relative to the workpath
params.RBM.continue = false;

params.RBM.inSamplesDims = [20 20 3]; % [x,y,#features]
params.RBM.inSamplesBorderBuffer = 0;

params.RBM.patchDim = 4; % after subsampling of input
params.RBM.hiddenSize = 100;
params.RBM.inputSubsampling = 2;

params.RBM.initWMax = 0.001;
params.RBM.initWGaussian = true;
params.RBM.initWScaleDistFcn = [];

params.RBM.saveinterval = 1;
params.RBM.saveintervalmat = 1;
params.RBM.saveintervalpng = 1;

params.RBM.topo = 0;
params.RBM.topoNghFcn = [];
params.RBM.topoNgh = [3 3];
params.RBM.topoGridDims = [20 20];
params.RBM.topoPeriodicBoundary = [true true];
params.RBM.topoEpsilon = 1e-2;

params.RBM.numberOfPatchReloads = 1000;
params.RBM.numberOfImagesPerPatchReload = 100;
params.RBM.numberOfPatchesPerReload = 500;
params.RBM.numberOfIterationsPerReload = 5;

params.RBM.validationSetsize = 0;
params.RBM.validationInterval = Inf;

params.RBM.lRate = {0.0001, 0.001, 0.01};			% LEARNING RATE
params.RBM.nGibbs = 1;				% # OF GIBBS SAMPLES (CONTRASTIVE DIVERGENCE)

params.RBM.sparsity = 0.01;		% TARGET HIDDEN UNIT SPARSITY
params.RBM.sparseGain = 0;			% GAIN ON THE LEARNING RATE FOR SPARSITY CONSTRAINTS
params.RBM.momentum = 0.9;			% (DEFAULT) GRADIENT MOMENTUM FOR WEIGHT UPDATES
params.RBM.wPenalty = 0;%.05;			% L2 WEIGHT PENALTY
params.RBM.beginAnneal = Inf;		% BEGIN SIMULUATED ANNEALING AFTER THIS # OF EPOCHS
params.RBM.beginWeightDecay = Inf;%1;	% BEGIN WEIGHT DECAY AFTER THIS # OF EPOCHS


paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


