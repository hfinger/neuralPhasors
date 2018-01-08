clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = '!(*thalia*|*erato*)';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
params.Gridjob.jobname = 'ConnectomeSimTest';
params.Gridjob.initRandStreamWithJobid = true;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '3';

params.ConnectomeSim.dataset = 4;
params.ConnectomeSim.subjId = [1:4 6:10 11:13 15 17:20];

params.ConnectomeSim.normRowBeforeHomotopic = 1;
params.ConnectomeSim.homotopic = num2cell(0:0.02:0.22);
params.ConnectomeSim.normRow = 1;

params.ConnectomeSim.model = 'SAR';
params.ConnectomeSim.normStd = true;
params.ConnectomeSim.k = num2cell(0.4:0.05:0.95);
params.ConnectomeSim.outFilenames = 'ConnectomeSimTest';
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);