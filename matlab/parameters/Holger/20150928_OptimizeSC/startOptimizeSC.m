clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'optim';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = [];%'nbp.q';
params.Gridjob.wc_host = '';

params.OptimizeSARSC.useFrob = false;
params.OptimizeSARSC.useRowRenorm = {'no','allowScaling'};
params.OptimizeSARSC.constrainPos = true;
params.OptimizeSARSC.k = 0.65;
params.OptimizeSARSC.lambda = {0, 0.1, 0.3, 1, 3, 10};
params.OptimizeSARSC.gamma = 1;
params.OptimizeSARSC.maxSteps = 262144;
params.OptimizeSARSC.savefolder = 'results';
params.OptimizeSARSC.unittest = false;
params.OptimizeSARSC.useRMSprop = true;
params.OptimizeSARSC.rmsPropDecay = 0.9;
params.OptimizeSARSC.initLearnRate = 1e-4;
params.OptimizeSARSC.saveAtIters = [0 2.^(0:18)];
params.OptimizeSARSC.breakAtCost = 0.999;
params.OptimizeSARSC.numRepeat = 10;
params.OptimizeSARSC.noiseAmount = 0:0.2:1;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);