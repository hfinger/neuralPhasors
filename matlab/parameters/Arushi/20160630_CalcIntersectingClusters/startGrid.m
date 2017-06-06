clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 15000;
params.Gridjob.jobname = 'CalcIntersectingClusters';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';
params.Gridjob.runOnlyJobIds =[];
params.Gridjob.wc_host = '!(*ramsauer*|*velma*|*fred*|*taygete*)'; %  | !shaggy';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope


params.CalcIntersectingClusters.split = num2cell(2:100:1000);
params.CalcIntersectingClusters.cosText = 'conn';
params.CalcIntersectingClusters.splitType = 'NonRec';
params.CalcIntersectingClusters.threshRange = [100];
params.CalcIntersectingClusters.subjRange = 1;
params.CalcIntersectingClusters.WeighFacRange = num2cell([0:0.1:1]);
            
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
