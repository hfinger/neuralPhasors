clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'Graclusjob2';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
% params.Gridjob.runOnlyJobIds = [6];
params.Gridjob.wc_host = '!(*ramsauer*|*daphne*|*orthosie*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.continue = true;

params.Graclusjob.clusterCount = num2cell([:40:1000]);
params.Graclusjob.WeighFactor = num2cell([0 0.2 0.4 0.5 0.6 0.8 1]);
params.Graclusjob.subjNum     = 1;
params.Graclusjob.decayParam    = -1;
params.Graclusjob.threshFactor  = 1;
params.Graclusjob.useCosineSim  = false;
params.Graclusjob.normBy        = 'sum';
params.Graclusjob.WholeNormText = 'WholeMax';
params.Graclusjob.CompSimPath   = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/';
params.Graclusjob.GraclusPath   =  '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_GraclusCut/';
            
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);