clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'GenerateClustConnmat';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';

params.GenerateClustConnmat.split = num2cell(2:1000);
params.GenerateClustConnmat.WeighingFactor = num2cell(0:0.1:0.9);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
