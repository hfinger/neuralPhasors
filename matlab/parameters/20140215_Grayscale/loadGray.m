clear paramsAll;

folderlist = { ...
'05june05_static_street_boston',...
'05june05_static_street_porter',...
'10feb04_static_cars_highland',...
'30may05_static_street_cambridge',...
'april21_static_outdoor_davis',...
'april21_static_outdoor_kendall',...
'barcelona_static_street',...
'boston_static_march',...
'dec_static_street',...
'madrid_static_street',...
'nov25_static_street',...
'nov6_static_outdoor',...
'nov7_static_outdoor',...
'oct6_static_outdoor',...
'paris_static_street',...
'static_barcelona_street_city_outdoor_2005',...
'static_harvard_outdoor_street',...
'static_outdoor_anchorage_alaska_usa'};
% Removed: 
% static_outdoor_city_laredo_spain
% static_dartmouth_hanover_june_2006
% static_nature_web_outdoor_animal
% static_newyork_city_urban
% static_outdoor_bay_area_submitted_alyosha_efros
% static_outdoor_bozeman_montana_usa
% static_outdoor_city_laredo_spain

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'labelMeInput';
params.LoadLabelMe.catName = folderlist;
params.LoadLabelMe.fileid = [];
params.LoadLabelMe.interpolateNearestNeighbor = true;
params.LoadLabelMe.convertToGray = true;
params.LoadLabelMe.outActFolder = 'labelMeInput';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'PreprocessZtrafo';
params.Gridjob.requiremf = 12000;
params.PreprocessZtrafo.inActFolder = 'labelMeInput';
params.PreprocessZtrafo.outWeightsFolder = 'zTrafoWeights';
params.PreprocessZtrafo.inNumChannels = 1;
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeZtransformed';
params.ApplyWeights.inActFolder = 'labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = [];
params.ApplyWeights.outActFolder = 'labelMeZtransformed';
params.ApplyWeights.weightsFile = 'zTrafoWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


