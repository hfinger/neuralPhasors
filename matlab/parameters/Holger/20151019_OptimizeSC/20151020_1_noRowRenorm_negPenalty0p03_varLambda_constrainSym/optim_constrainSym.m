clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'optim_constrainSym';
params.Gridjob.initRandStreamWithSeed = 0;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = [];%'nbp.q';
params.Gridjob.wc_host = '';

params.OptimizeSARSC.constrainSym = true;
params.OptimizeSARSC.useFrob = false;
params.OptimizeSARSC.useRowRenorm = 'no';
params.OptimizeSARSC.constrainPos = false;
params.OptimizeSARSC.k = 0.65;
params.OptimizeSARSC.lambda = {0, 0.001, 0.003, 0.01};
params.OptimizeSARSC.gamma = 1;
params.OptimizeSARSC.negPenalty = 0.03;
params.OptimizeSARSC.maxSteps = 262144;
params.OptimizeSARSC.savefolder = 'results';
params.OptimizeSARSC.unittest = false;
params.OptimizeSARSC.useRMSprop = true;
params.OptimizeSARSC.rmsPropDecay = 0.9;
params.OptimizeSARSC.initLearnRate = 1e-4;
params.OptimizeSARSC.saveAtIters = [0 2.^(0:18)];
params.OptimizeSARSC.breakAtCost = 0.999999;
params.OptimizeSARSC.noiseAmount = num2cell(0:0.25:1);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);