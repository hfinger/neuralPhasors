clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'ProbtrackXallsubjects';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';
params.Gridjob.queue = 'nbp.q';

params.ProbtrackXallsubjects.split = num2cell(50);
params.ProbtrackXallsubjects.numSubjects = 22;
params.ProbtrackXallsubjects.splitPerSubject = 50;


paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
