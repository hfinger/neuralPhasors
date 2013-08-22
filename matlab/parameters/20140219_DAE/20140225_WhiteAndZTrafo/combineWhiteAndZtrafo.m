clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.wc_host = [];
params.Gridjob.jobname = 'concat';
params.Gridjob.initRandStreamWithSeed = 12345;
params.ActConcat.inActFolder1 = '../../20131217_ReLuDAE/20140204_InputZtrafo/labelMeZtransformed';
params.ActConcat.inActFilenames1 = 'act.*.mat';
params.ActConcat.inActFolder2 = '../../20131217_ReLuDAE/20131220_MoreImages/labelMeWhite';
params.ActConcat.inActFilenames2 = 'act.*.mat';
params.ActConcat.outActFolder = 'concatInputs'; %relative to the workpath
paramsAll{1} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);