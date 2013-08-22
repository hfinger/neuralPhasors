
%%
clear params;

params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'calcInputCov';

params.FeatureCovariance.inActFolder = '../20130318_AE400/labelMeInput';
params.FeatureCovariance.inActFilenames = 'act.*.mat';
params.FeatureCovariance.fileid = [];
params.FeatureCovariance.outCovFolder = 'InputCorrCov'; %relative to the workpath
params.FeatureCovariance.maxCovLengthDim1 = 149;
params.FeatureCovariance.maxCovLengthDim2 = 199;
params.FeatureCovariance.numSamplesPerImage = 1;
params.FeatureCovariance.borderBuffer = 0; %exclude these pixels around the image
params.FeatureCovariance.saveCorr = false;
params.FeatureCovariance.saveCov = false;
params.FeatureCovariance.saveVar = true;

gridjobs = Gridjob(params);
start(gridjobs);
