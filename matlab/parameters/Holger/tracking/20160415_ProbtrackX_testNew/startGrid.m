clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2500;
params.Gridjob.jobname = 'ProbtrackXsubj3and20only';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.runOnlyJobIds = [];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*kalyke*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;

params.ProbtrackXallsubjects.split = num2cell(1:100);
params.ProbtrackXallsubjects.splitPerSubject = 100;
params.ProbtrackXallsubjects.subjectId = num2cell([3 20]);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);