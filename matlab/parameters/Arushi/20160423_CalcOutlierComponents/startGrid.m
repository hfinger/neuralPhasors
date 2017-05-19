clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 1000;
params.Gridjob.jobname = 'CalcOutlierComponents';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';

params.CalcOutlierComponents.split = num2cell(287);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
