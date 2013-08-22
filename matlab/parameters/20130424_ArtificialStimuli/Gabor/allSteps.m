clear paramsAll;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.jobname = 'imageInput';
params.Gridjob.wc_host = '';
params.LoadImage.filepath = {...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3pt.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3ptBroken.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3ptColorGradientOnGray.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/continousLines3ptColorGradient.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/cube.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/Kanizsa_triangle.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/Kanizsa_triangle_only_dots.png',...
  '/net/store/nbp/phasesim/databases/ArtificialStimuli/Kanizsa_triangle_only_triangle.png'};
params.LoadImage.outActFolder = 'imageInput';
paramsAll{1} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 13000;
params.Gridjob.jobname = 'imageWhite';
params.Gridjob.wc_host = '';
params.ApplyWeights.inActFolder = 'imageInput';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = num2cell(1:9);
params.ApplyWeights.outActFolder = 'imageWhite';
params.ApplyWeights.weightsFile = '../../20130409_GaborNoWhiteningWithSigmoidWithZtrafo5timesWithSubsampling/labelMeWhiteningWeights/weights.mat';
params.ApplyWeights.convType = 'same';
paramsAll{2} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActNotRectified';
params.Gridjob.wc_host = '';
params.ApplyWeights.inActFolder = 'imageWhite';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:9;
params.ApplyWeights.outActFolder = 'layer1ActNotRectified';
params.ApplyWeights.weightsFile = '../../20130409_GaborNoWhiteningWithSigmoidWithZtrafo5timesWithSubsampling/labelMeAE100weights/forwConn.mat';
params.ApplyWeights.convType = 'same';
paramsAll{3} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 12000;
params.Gridjob.jobname = 'layer1ActRectified';
params.Gridjob.wc_host = '';
params.ApplyWeights.inActFolder = 'layer1ActNotRectified';
params.ApplyWeights.inActFilenames = 'act.*.mat';
params.ApplyWeights.fileid = 1:9;
params.ApplyWeights.outActFolder = 'layer1ActRectified';
params.ApplyWeights.weightsFile = [];
params.ApplyWeights.actFcn2 = @(x) feval(@(x2) bsxfun(@rdivide,x2,0.001+sum(x2,3)), max(0,bsxfun(@minus,x,mean(x,3))) ); % subtract mean and rectify and then normalize to sum=1
params.ApplyWeights.convType = 'same';
paramsAll{4} = params;

clear params;
params.Gridjob.runLocal = false;
params.Gridjob.requiremf = 4000;
params.Gridjob.jobname = 'layer1Phase';
params.Gridjob.wc_host = '';
params.PhaseSimulation.inActFolder = 'layer1ActRectified';
params.PhaseSimulation.inActFilenames = 'act.*.mat';
params.PhaseSimulation.inFileid = num2cell(1:9);
params.PhaseSimulation.inCellid = 1;
params.PhaseSimulation.inConnFilename = '../../20130409_GaborNoWhiteningWithSigmoidWithZtrafo5timesWithSubsampling/layer1Conn/weights.mat';
params.PhaseSimulation.inPhaseFilename = [];
params.PhaseSimulation.outPhaseFolder = 'layer1Phase';
params.PhaseSimulation.noiseLevel = 0;
params.PhaseSimulation.noiseEMAconst = 0;
params.PhaseSimulation.tmax = 50;
params.PhaseSimulation.dt = 1;
params.PhaseSimulation.fixedPhaseDelay = 0;
params.PhaseSimulation.odeSolver = 'ode4';
params.PhaseSimulation.weightAll = 3;
params.PhaseSimulation.weightInh = 1;
params.PhaseSimulation.weightExc = 1;
params.PhaseSimulation.saveintervalPhase = 5;
params.PhaseSimulation.saveintervalMeanPhase = 5;
params.PhaseSimulation.saveintervalMeanWeightedPhase = 1;
params.PhaseSimulation.plotPhase = false;
params.PhaseSimulation.maxdphase = 0.5;
paramsAll{5} = params;

clear params;
gridjobs = Gridjob(paramsAll);
start(gridjobs);


