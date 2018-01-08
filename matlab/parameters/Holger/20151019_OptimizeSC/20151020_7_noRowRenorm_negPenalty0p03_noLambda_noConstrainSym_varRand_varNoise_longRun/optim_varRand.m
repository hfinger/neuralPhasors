clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'optim_varRand';
params.Gridjob.initRandStreamWithSeed = {0, 1};
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = [];%'nbp.q';
params.Gridjob.wc_host = '';

params.OptimizeSARSC.constrainSym = false;
params.OptimizeSARSC.useFrob = false;
params.OptimizeSARSC.useRowRenorm = 'no';
params.OptimizeSARSC.constrainPos = false;
params.OptimizeSARSC.k = 0.65;
params.OptimizeSARSC.lambda = 0;
params.OptimizeSARSC.gamma = 1;
params.OptimizeSARSC.negPenalty = 0.03;
params.OptimizeSARSC.maxSteps = 2^20;
params.OptimizeSARSC.savefolder = 'results';
params.OptimizeSARSC.unittest = false;
params.OptimizeSARSC.useRMSprop = true;
params.OptimizeSARSC.rmsPropDecay = 0.9;
params.OptimizeSARSC.initLearnRate = 1e-4;
params.OptimizeSARSC.saveAtIters = [0 2.^(0:20)];
params.OptimizeSARSC.breakAtCost = 0.99999999;
params.OptimizeSARSC.noiseAmount = num2cell(0:0.1:1);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);