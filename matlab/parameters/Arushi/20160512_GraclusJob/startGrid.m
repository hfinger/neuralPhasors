clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2500;
params.Gridjob.jobname = 'Graclusjob2';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
% params.Gridjob.runOnlyJobIds = [39];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*kalyke*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;

params.Graclusjob.WeighFactor = num2cell(1);
params.Graclusjob.clusterCount = num2cell(2:30:1000);
params.Graclusjob.subjNum     = 1;
params.Graclusjob.decayParam    = -1;
params.Graclusjob.threshFactor  = 100;
params.Graclusjob.useCosineSim  = false;
params.Graclusjob.normBy        = 'sum';
params.Graclusjob.WholeNormText = 'WholeMax';
params.Graclusjob.CompSimPath   = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/';
params.Graclusjob.GraclusPath   =  '/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/';
            
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);