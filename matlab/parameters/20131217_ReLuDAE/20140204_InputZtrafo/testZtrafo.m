clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'PreprocessZtrafoTest';
params.Gridjob.requiremf = 12000;
params.PreprocessZtrafo.inActFolder = 'labelMeZtransformed';
params.PreprocessZtrafo.outWeightsFolder = 'zTrafoWeightsTest';
params.PreprocessZtrafo.inNumChannels = 3;
paramsAll{1} = params;


clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


