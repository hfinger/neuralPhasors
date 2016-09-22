clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2500;
params.Gridjob.jobname = 'ProbtrackXallsubjects';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.runOnlyJobIds = [520, 521, 524, 528, 529, 530, 543];
params.Gridjob.wc_host = '!(*ramsauer*) |!(*daphne*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;



params.ProbtrackXallsubjects.split = num2cell(1:1100);
params.ProbtrackXallsubjects.numSubjects = 22;
params.ProbtrackXallsubjects.splitPerSubject = 50;


paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);