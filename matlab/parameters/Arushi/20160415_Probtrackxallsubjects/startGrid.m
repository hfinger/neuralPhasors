clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 2500;
params.Gridjob.jobname = 'ProbtrackXallsubjects';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '2';
params.Gridjob.queue = 'nbp.q';
% params.Gridjob.runOnlyJobIds = [520, 521, 524, 528, 529, 530, 543];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*kalyke*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;

params.ProbtrackXallsubjects.split = num2cell(1:100);
params.ProbtrackXallsubjects.splitPerSubject = 100;
params.ProbtrackXallsubjects.subjectId = num2cell(1);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);