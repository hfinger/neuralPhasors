clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'AE100SegEval';
params.Gridjob.combParallel = true;
params.SegmentationEval.inActFolder = 'AE100layer1Act';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'AE100layer1PhaseAll';
params.SegmentationEval.outSegEvalFolder = 'AE100Layer1SegEval';
params.SegmentationEval.catName = '05june05_static_street_boston';
params.SegmentationEval.fileid = 1:50;
params.SegmentationEval.minPixelPerSeg = 12.5;
params.SegmentationEval.maxPixelPerSeg = Inf;
params.SegmentationEval.borderSize=0;
params.SegmentationEval.numPairs=10000; % set numPairs=0 to calc only probSegSync
params.SegmentationEval.numSubpop=100;
params.SegmentationEval.numSamplesInSubpop=1000;
params.SegmentationEval.phasePerPixel=false;
params.SegmentationEval.time={0,5,10,15,20,25,30};
params.SegmentationEval.initRandstream = true;
params.SegmentationEval.pairOnlyWithOneNonmatching = true;
params.SegmentationEval.verboseLevel = 0;
params.SegmentationEval.simOtherVarParam = [];
params.SegmentationEval.skipMissingImg = false;
params.SegmentationEval.borderEvalNumPerImg = 100;
params.SegmentationEval.borderEvalNghSize = 10;
params.SegmentationEval.borderAngleSigma = 3;
params.SegmentationEval.borderAngleUseActivityWeighting = false;
params.SegmentationEval.borderSyncUseActivityWeighting = true;
params.SegmentationEval.useHueInsteadPhase = false;
params.SegmentationEval.resizeToX = 200;
params.SegmentationEval.resizeToY = 150;
params.SegmentationEval.cropX = [];
params.SegmentationEval.cropY = [];
params.SegmentationEval.doCalcSegSync = true;
params.SegmentationEval.doCalcBorderSync = true;
params.SegmentationEval.doCalcBorderAngle = true;

gridjobs = Gridjob(params);
start(gridjobs);
