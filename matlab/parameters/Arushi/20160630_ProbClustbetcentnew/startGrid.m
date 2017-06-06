clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 10000;
params.Gridjob.jobname = 'ProbClustBetCentnew';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = true;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.runOnlyJobIds =[];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*velma*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope


params.ProbClustBetCentnew.split = num2cell([2:100:1000]);
params.ProbClustBetCentnew.cosText = 'conn'; %cos
params.ProbClustBetCentnew.splitType = 'NonRec'; %'NonRec'
params.ProbClustBetCentnew.threshRange =  [1];
params.ProbClustBetCentnew.WeighFacRange = num2cell([0.5]);

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
