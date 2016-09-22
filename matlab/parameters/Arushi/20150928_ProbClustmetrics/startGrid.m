clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 1000;
params.Gridjob.jobname = 'ProbClustbetcent';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = 'nbp.q';

params.ProbClustmetrics.split = num2cell(1:1276);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
