clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 10000;
params.Gridjob.jobname = 'ProbClustBetCentnew';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*velma*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope


params.ProbClustBetCentnew.split = num2cell(0:0.1:1);
params.ProbClustBetCentnew.cosText = 'conn'; %cos
params.ProbClustBetCentnew.splitType = 'Rec'; %'NonRec'
params.ProbClustBetCentnew.threshRange =  [100,10,5,1];

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
