clear paramsAll;

clear params;
params.Gridjob.runLocal = true;
params.Gridjob.jobname = 'labelMeInput';
params.LoadImage.resizeToX = 400;
params.LoadImage.resizeToY = 300;
params.LoadImage.doResize = true;
params.LoadImage.filepath = cellfun(@(x) num2str(x, '/net/store/nbp/phasesim/databases/PhaseFeat/labelme/image_%02d.jpg'), num2cell(1:100),'UniformOutput',false);
params.LoadImage.outActFolder = 'labelMeInput'; %relative to the workpath      
params.LoadImage.useFilenameAsDir = true;
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'labelMeWhite';
params.ApplyWeights.inActFolder = 'labelMeInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:100;
params.ApplyWeights.outActFolder = 'labelMeWhite';
params.ApplyWeights.weightsFile = '../../20130726_Paper/Autoencoder/labelMeWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActNotRectified';
params.ApplyWeights.inActFolder = 'labelMeWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:100;
params.ApplyWeights.outActFolder = 'layer1ActNotRectified';
params.ApplyWeights.weightsFile = '../../20130726_Paper/Autoencoder/labelMeAE100weights/forwConn.mat';
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActRectified';
params.ApplyWeights.inActFolder = 'layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:100;
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 2000;
params.Gridjob.jobname = 'phaseSim';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:100);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = '../../20130726_Paper/Autoencoder/layer1ConnManySamplesFDR0p05maxdx32/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'phaseSim';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 30;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 10;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 5;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 2*pi/4;
paramsAll{5} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1PhaseFeat';
params.PhaseFeat.inActFolder = 'layer1ActRectified';
params.PhaseFeat.inActFilenames = 'act.*.mat';
params.PhaseFeat.inPhaseFolder = 'phaseSim';
params.PhaseFeat.inPhaseFilenames = 'phaseIter30.mat';
params.PhaseFeat.fileid = 1:100;
params.PhaseFeat.outPhaseFeatFolder = 'layer1PhaseFeat';
params.PhaseFeat.filtRadius = {1,2,3,5};
params.PhaseFeat.coherenceFactor = 1;
paramsAll{6} = params;

clear params;
allfixmat = load('/net/store/nbp/PhaseFeat/labelme/allfixmat.mat');
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'layer1FixationEval';
params.EvaluateFixations.inPhaseFeatFolder = {'layer1PhaseFeat/1','layer1PhaseFeat/2','layer1PhaseFeat/3','layer1PhaseFeat/4'};
params.EvaluateFixations.inPhaseFeatFilenames = 'act.*.mat';
params.EvaluateFixations.fileid = 1:100;
params.EvaluateFixations.inFixationMat = allfixmat.fixmat;
params.EvaluateFixations.inFixationFileid = 1:100;
params.EvaluateFixations.inFixationScaleX = 1/6.4;
params.EvaluateFixations.inFixationScaleY = 1/6.4;
params.EvaluateFixations.inFixationOffsetX = 320;
params.EvaluateFixations.inFixationOffsetY = 0;
params.EvaluateFixations.outFolder = 'layer1FixationEval';
paramsAll{7} = params;

clear params;
allfixmat = load('/net/store/nbp/PhaseFeat/labelme/allfixmat.mat');
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 5000;
params.Gridjob.jobname = 'layer1FixationEvalNonmatch';
params.EvaluateFixations.inPhaseFeatFolder = {'layer1PhaseFeat/1','layer1PhaseFeat/2','layer1PhaseFeat/3','layer1PhaseFeat/4'};
params.EvaluateFixations.inPhaseFeatFilenames = 'act.*.mat';
params.EvaluateFixations.fileid = 1:100;
params.EvaluateFixations.inFixationMat = allfixmat.fixmat;
params.EvaluateFixations.inFixationFileid = [2:100 1];
params.EvaluateFixations.inFixationScaleX = 1/6.4;
params.EvaluateFixations.inFixationScaleY = 1/6.4;
params.EvaluateFixations.inFixationOffsetX = 320;
params.EvaluateFixations.inFixationOffsetY = 0;
params.EvaluateFixations.outFolder = 'layer1FixationEvalNonmatch';
paramsAll{8} = params;

clear params;
gridjobs = Gridjob(paramsAll(7:8));
start(gridjobs);


