clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 27000;
params.Gridjob.jobname = 'ProbtrackXappendnorm';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = 'nbp.q';

params.ProbtrackXappendnorm.split = num2cell(1:1);
params.ProbtrackXappendnorm.numberPerSplit = 1000;
params.ProbtrackXappendnorm.numberRegions = 52442;


paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
