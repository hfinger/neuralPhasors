clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'lMLayer1SegEvalNewManySamplesHueBorder0AllImg';
params.Gridjob.combParallel = true;
params.SegmentationEval.inActFolder = 'labelMeActLayer1';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'labelMeLayer1PhaseAll';
params.SegmentationEval.outSegEvalFolder = 'labelMeLayer1SegEvalNewManySamplesHueBorder0AllImg';
params.SegmentationEval.catName = '05june05_static_street_boston';
params.SegmentationEval.fileid = 1:185;
params.SegmentationEval.minPixelPerSeg = 50;
params.SegmentationEval.maxPixelPerSeg = Inf;
params.SegmentationEval.borderSize=0;
params.SegmentationEval.numPairs=[]; % set numPairs=0 to calc only probSegSync
params.SegmentationEval.numSubpop=1000;
params.SegmentationEval.numSamplesInSubpop=25;
params.SegmentationEval.phasePerPixel=false;
params.SegmentationEval.time=0;
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
params.SegmentationEval.useHueInsteadPhase = true;
params.SegmentationEval.resizeToX = 400;
params.SegmentationEval.resizeToY = 300;
params.SegmentationEval.cropX = [];
params.SegmentationEval.cropY = [];
params.SegmentationEval.doCalcSegSync = true;
params.SegmentationEval.doCalcBorderSync = true;
params.SegmentationEval.doCalcBorderAngle = true;

gridjobs = Gridjob(params);
start(gridjobs);
