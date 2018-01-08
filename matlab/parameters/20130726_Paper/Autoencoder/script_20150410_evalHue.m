
clear params;
params.Gridjob.runLocal = true;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'segEvalHue';
params.Gridjob.combParallel = true;
params.SegmentationEval.inActFolder = 'layer1ActRectified';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'img50PhaseManySamplesFDR0p05maxdx32SmallDt';
params.SegmentationEval.outSegEvalFolder = 'segEvalHue';
params.SegmentationEval.catName = '05june05_static_street_boston';
params.SegmentationEval.fileid = [1:50];
params.SegmentationEval.segSizeBinEdges = [exp(0:11)];
params.SegmentationEval.minPixelPerSeg = 1;
params.SegmentationEval.maxPixelPerSeg = Inf;
params.SegmentationEval.borderSize=0;
params.SegmentationEval.numPairs=10000;
params.SegmentationEval.numSubpop=100;
params.SegmentationEval.numSamplesInSubpop=1000;
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
params.SegmentationEval.resizeToX = 200;
params.SegmentationEval.resizeToY = 150;
params.SegmentationEval.cropX = [];
params.SegmentationEval.cropY = [];
params.SegmentationEval.doCalcSegSync = true;
params.SegmentationEval.doCalcBorderSync = true;
params.SegmentationEval.doCalcBorderAngle = true;
paramsAll{1} = params;


clear params;
gridjobs = Gridjob(paramsAll(1));
start(gridjobs);

