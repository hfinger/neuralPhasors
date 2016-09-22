clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 1000;
params.Gridjob.jobname = 'ProbClustbetcent';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = 'nbp.q';

params.ProbClustbetcentsub18.split = num2cell(1013:1044);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
