clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 6000;
params.Gridjob.jobname = 'ProbtrackX';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = 'nbp.q';

params.ProbtrackX.split = num2cell(1:53);
params.ProbtrackX.numberPerSplit = 1000;
params.ProbtrackX.numberRegions = 52442;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);