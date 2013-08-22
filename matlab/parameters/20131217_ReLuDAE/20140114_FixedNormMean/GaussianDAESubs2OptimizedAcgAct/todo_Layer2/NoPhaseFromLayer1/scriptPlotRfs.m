clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1PlotRfImagePatches';
params.PlotRfImagePatches.inActImageFolder = '../../20130726_Paper/Autoencoder/labelMeInput';
params.PlotRfImagePatches.inActImageFilenames = 'act.*.mat';
params.PlotRfImagePatches.inActFolder = 'layer1ActRectified';
params.PlotRfImagePatches.inActFilenames = 'act.*.mat';
params.PlotRfImagePatches.plotNumPatches = 10;
params.PlotRfImagePatches.plotPatchHalfLength = 6; %measured in units of the image layer
params.PlotRfImagePatches.excludePatchHalfLength = 6; %measured in units of the activity layer
params.PlotRfImagePatches.outPatchesFolder = 'layer1PlotRfImagePatches';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer2PlotRfImagePatches';
params.PlotRfImagePatches.inActImageFolder = '../../20130726_Paper/Autoencoder/labelMeInput';
params.PlotRfImagePatches.inActImageFilenames = 'act.*.mat';
params.PlotRfImagePatches.inActFolder = 'layer2ActRectified';
params.PlotRfImagePatches.inActFilenames = 'act.*.mat';
params.PlotRfImagePatches.plotNumPatches = 10;
params.PlotRfImagePatches.plotPatchHalfLength = 12; %measured in units of the image layer
params.PlotRfImagePatches.excludePatchHalfLength = 6; %measured in units of the activity layer
params.PlotRfImagePatches.outPatchesFolder = 'layer2PlotRfImagePatches';
paramsAll{2} = params;

clear params;
gridjobs = Gridjob(paramsAll(1:2));
start(gridjobs);


