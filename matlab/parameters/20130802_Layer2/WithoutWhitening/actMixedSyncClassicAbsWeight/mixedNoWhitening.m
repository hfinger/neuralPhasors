clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'initL2';
params.ApplyWeightsSyncAndClassic.inActFolder = '../../NoPhaseFromLayer1/layer1ActRectified';
params.ApplyWeightsSyncAndClassic.inActFilenames = 'act.*.mat';
params.ApplyWeightsSyncAndClassic.inPhaseFolder = '../../../20130726_Paper/Autoencoder/img50PhaseManySamplesFDR0p05maxdx32SmallDt';
params.ApplyWeightsSyncAndClassic.inPhaseFilenames = 'phaseIter30.mat';
params.ApplyWeightsSyncAndClassic.fileid = 1:50;
params.ApplyWeightsSyncAndClassic.outActFolder = 'initActL2';
params.ApplyWeightsSyncAndClassic.outPhaseFolder = 'initPhaseL2';
params.ApplyWeightsSyncAndClassic.weightsFile = '../temp_layer2AE100weights/forwConnIter120.mat';
params.ApplyWeightsSyncAndClassic.weightSyncTerm = {0, 0.5, 1};
params.ApplyWeightsSyncAndClassic.convType = 'same';
params.ApplyWeightsSyncAndClassic.useAbsWeight = true;
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);
