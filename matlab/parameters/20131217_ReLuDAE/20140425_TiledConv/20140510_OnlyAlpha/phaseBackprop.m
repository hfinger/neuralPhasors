clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'phaseBackpropL1Input';

params.ApplyWeightsSyncAndClassic.inActFolder = 'layer1ActNotRectified/p1010736.jpg'; %relative to the workpath
params.ApplyWeightsSyncAndClassic.inActFilenames = 'act1.mat';
params.ApplyWeightsSyncAndClassic.inPhaseFolder = 'phase/2/p1010736.jpg'; %relative to the workpath
params.ApplyWeightsSyncAndClassic.inPhaseFilenames = 'phaseIter30.mat';
params.ApplyWeightsSyncAndClassic.fileid = [];
params.ApplyWeightsSyncAndClassic.outActFolder = []; %relative to the workpath
params.ApplyWeightsSyncAndClassic.outPhaseFolder = 'phaseBackpropL1Input'; %relative to the workpath
params.ApplyWeightsSyncAndClassic.outStatsFolder = []; %relative to the workpath
params.ApplyWeightsSyncAndClassic.weightsFile = 'ReLuDAE/backConn.mat'; %relative to the workpath
params.ApplyWeightsSyncAndClassic.plotPdf = false; %if true then plot a pdf
params.ApplyWeightsSyncAndClassic.weightSyncTerm = 0.5; % between 0 (classic) and 1 (sync) to weight between classic and sync term
params.ApplyWeightsSyncAndClassic.useAbsWeight = true; % TODO check if this makes sense
params.ApplyWeightsSyncAndClassic.calcStatsBins = 0;
% If the following parameters are set, they will overwrite the specifications in the weightsFile
params.ApplyWeightsSyncAndClassic.inputSubsampling = [];
params.ApplyWeightsSyncAndClassic.shiftOutputdims = true;
params.ApplyWeightsSyncAndClassic.outputUpSampling = 2;
params.ApplyWeightsSyncAndClassic.convType = 'same'; %{valid},same,full
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


