clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 10000;
params.Gridjob.jobname = 'GenerateClustConnmat';
params.Gridjob.initRandStreamWithJobid = false;
params.Gridjob.continue = false;
params.Gridjob.requiredThreads = '1';
params.Gridjob.queue = 'nbp.q';

params.GenerateClustConnmat.WeighFactor = num2cell(0:0.1:1);
params.GenerateClustConnmat.clusterCount = (2:1000);
params.GenerateClustConnmat.subjNum     = 1;
params.GenerateClustConnmat.decayParam    = -1;
params.GenerateClustConnmat.useCosineSim  = false;
params.GenerateClustConnmat.recursiveSplit = false;
params.GenerateClustConnmat.normBy        = 'sum';
params.GenerateClustConnmat.WholeNormText = 'WholeMax';
params.GenerateClustConnmat.CompSimPath   = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/';
params.GenerateClustConnmat.GraclusPath   =  '/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/';
params.GenerateClustConnmat.PostProcessPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/';
params.GenerateClustConnmat.threshRange   = [100,10,5,1];

paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
