clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 1000;
params.Gridjob.jobname = 'ProbClustbetcent';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = 'nbp.q';

params.ProbClustbetcent.split = num2cell(693);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
