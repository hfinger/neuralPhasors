clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 1000;
params.Gridjob.jobname = 'CalcOutlierComponents';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';

params.CalcOutlierComponents.WeighingFactor = num2cell([0:0.1:0.4 0.6:0.1:1]);
params.CalcOutlierComponents.recursiveSplit = 'NonRec';
params.CalcOutlierComponents.cosText        = 'conn';
params.CalcOutlierComponents.clustRange     = 2:1000;
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
