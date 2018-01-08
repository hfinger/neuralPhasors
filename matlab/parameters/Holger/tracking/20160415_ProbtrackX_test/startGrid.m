clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 2500;
params.Gridjob.jobname = 'ProbtrackXsubj01only';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.runOnlyJobIds = [];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;

params.ProbtrackXallsubjects.split = num2cell(1:20);
params.ProbtrackXallsubjects.numSubjects = 1;
params.ProbtrackXallsubjects.splitPerSubject = 20;

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);