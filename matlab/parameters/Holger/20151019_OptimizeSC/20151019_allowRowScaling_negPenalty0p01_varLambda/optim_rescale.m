clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'optim_rescale';
params.Gridjob.initRandStreamWithSeed = 0;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = [];%'nbp.q';
params.Gridjob.wc_host = '';

params.OptimizeSARSC.useFrob = false;
params.OptimizeSARSC.useRowRenorm = 'allowScaling';
params.OptimizeSARSC.constrainPos = false;
params.OptimizeSARSC.k = 0.65;
params.OptimizeSARSC.lambda = {0, 0.01, 0.03, 0.1, 0.3, 1};
params.OptimizeSARSC.gamma = 1;
params.OptimizeSARSC.negPenalty = 0.01;
params.OptimizeSARSC.maxSteps = 262144;
params.OptimizeSARSC.savefolder = 'results';
params.OptimizeSARSC.unittest = false;
params.OptimizeSARSC.useRMSprop = true;
params.OptimizeSARSC.rmsPropDecay = 0.9;
params.OptimizeSARSC.initLearnRate = 1e-4;
params.OptimizeSARSC.saveAtIters = [0 2.^(0:18)];
params.OptimizeSARSC.breakAtCost = 0.999999;
params.OptimizeSARSC.numRepeat = 1;
params.OptimizeSARSC.noiseAmount = num2cell(0:0.25:1);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);