clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 10000;
params.Gridjob.jobname = 'ProbClustmetricsnew';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*velma*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
% params.Gridjob.runOnlyJobIds = [2];


params.ProbClustmetricsnew.split = num2cell(0.5);
params.ProbClustmetricsnew.cosText = 'cos'; %cos
params.ProbClustmetricsnew.splitType = 'Rec'; %'NonRec'
params.ProbClustmetricsnew.threshRange =  [1];

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
