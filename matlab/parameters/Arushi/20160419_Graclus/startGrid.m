clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2500;
params.Gridjob.jobname = 'Graclusjob1';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
% params.Gridjob.runOnlyJobIds = [520, 521, 524, 528, 529, 530, 543];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*kalyke*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;

params.Graclusjob.WeighFactor = num2cell(0:0.1:0.9);
params.Graclusjob.clusterCount = num2cell(201:20:500);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);