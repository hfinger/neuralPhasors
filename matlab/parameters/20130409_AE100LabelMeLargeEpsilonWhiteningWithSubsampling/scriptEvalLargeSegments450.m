clear paramsAll;

% clear params;
% params.Gridjob.runLocal = false;
% params.Gridjob.requiremf = 3000;
% params.Gridjob.jobname = 'layer1PhaseAll';
% params.PhaseSimulation.inActFolder = 'layer1ActRectified';
% params.PhaseSimulation.inActFilenames = 'act.*.mat';
% params.PhaseSimulation.inFileid = num2cell(1:185);
% params.PhaseSimulation.inCellid = 1;
% params.PhaseSimulation.inConnFilename = 'layer1Conn/weights.mat';
% params.PhaseSimulation.inPhaseFilename = [];
% params.PhaseSimulation.outPhaseFolder = 'layer1PhaseAll';
% params.PhaseSimulation.noiseLevel = 0;
% params.PhaseSimulation.noiseEMAconst = 0;
% params.PhaseSimulation.tmax = 30;
% params.PhaseSimulation.dt = 1;
% params.PhaseSimulation.fixedPhaseDelay = 0;
% params.PhaseSimulation.odeSolver = 'ode4';
% params.PhaseSimulation.weightAll = 30;
% params.PhaseSimulation.weightInh = 1;
% params.PhaseSimulation.weightExc = 1;
% params.PhaseSimulation.saveintervalPhase = 5;
% params.PhaseSimulation.saveintervalMeanPhase = 5;
% params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
% params.PhaseSimulation.plotPhase = false;
% params.PhaseSimulation.maxdphase = 0.5;
% paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 3000;
params.Gridjob.jobname = 'layer1EvalLargeSeg450';
params.Gridjob.combParallel = true;
params.SegmentationEval.inActFolder = 'layer1ActRectified';
params.SegmentationEval.inActFilenames = 'act1.mat';
params.SegmentationEval.inPhaseFolder = 'layer1PhaseAll';
params.SegmentationEval.outSegEvalFolder = 'layer1EvalLargeSeg450';
params.SegmentationEval.catName = '05june05_static_street_boston';
params.SegmentationEval.fileid = [1:115 117:147 149:185];
params.SegmentationEval.minPixelPerSeg = 450;
params.SegmentationEval.maxPixelPerSeg = Inf;
params.SegmentationEval.borderSize=0;
params.SegmentationEval.numPairs=10000;
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
paramsAll{2} = params;

gridjobs = Gridjob(paramsAll{2});
start(gridjobs);
