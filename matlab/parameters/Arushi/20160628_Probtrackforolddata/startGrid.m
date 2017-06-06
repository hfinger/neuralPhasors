clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2500;
params.Gridjob.jobname = 'Probtrackforolddata';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '2';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.runOnlyJobIds = [8];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*kalyke*|*velma*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;

params.Probtrackforolddata.split = num2cell([1:100]);
params.Probtrackforolddata.splitPerSubject = 100;
params.Probtrackforolddata.subjectId = num2cell(1);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);