classdef Autoencoder < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
    tempsavepath
    forwardRfSize
    backwardRfSize
    dimW1
    dimW2
    nghMat
    randstreamStateBeforeEval
    randstreamStateBeforeIter
    randstreamStateBeforeReload
    randstreamStateBeforeValidation
    
    tiledConvF
    tiledConvB
    tiledConvBInv
    
    theta
    nextReloadId = 0
    gradStates = []
    inputApplyConn = []
    EMA_z2pos
    EMA_y2_pairs
    EMA_y2_square
    EMA_y2
    currentNumberOfUnits
    
    log = struct();
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = Autoencoder(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.Autoencoder.inActFolder = 'whitenedData'; %relative to the workpath
      this.params.Autoencoder.inActFilenames = 'act.*.mat';
      this.params.Autoencoder.inputApplyConn = [];
      this.params.Autoencoder.outWeightsFolder = 'AEWeights'; %relative to the workpath
      this.params.Autoencoder.continue = false;
      this.params.Autoencoder.continueBatch = false; % !!!!!!! not working correctly! TODO: reorder dimensions of W!
      this.params.Autoencoder.continueBatchInWeightsFolder = 'AEWeights';
      this.params.Autoencoder.continueBatchInBackConnFilenames = 'backConnIter50.mat';
      this.params.Autoencoder.continueBatchInForwConnFilenames = 'forwConnIter50.mat';
      this.params.Autoencoder.loadPatchesAndTheta = [];

      % INIT PARAMS:
      this.params.Autoencoder.inSamplesDims = [32 32 3 500]; % [x,y,#features,#samples]
      this.params.Autoencoder.inSamplesBorderBuffer = 0;
      this.params.Autoencoder.patchDimForward = 8;
      this.params.Autoencoder.patchDimBackward = 8;
      this.params.Autoencoder.hiddenSize = 400;
      this.params.Autoencoder.useTiledConv = false; %if true: tileSize=patchDimForward-1, fOut=hiddenSize      
      this.params.Autoencoder.inputSubsampling = 1;
      this.params.Autoencoder.forwInitMaxWeight = 1;
      this.params.Autoencoder.forwInitScaleDistFcn = [];
      this.params.Autoencoder.backInitMaxWeight = 1;
      this.params.Autoencoder.backInitScaleDistFcn = [];
      this.params.Autoencoder.initWGaussian = false;
      this.params.Autoencoder.initWExponent = 1;
      
      % PROCESSING STAGES:
      this.params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
      this.params.Autoencoder.maskingNoiseFractionRescale = true;        % if inputs are rescaled to the variance before masking
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      this.params.Autoencoder.gaussianNoiseCorrelated = false;
      this.params.Autoencoder.gaussianNoiseCorrelatedUniform = false;

      this.params.Autoencoder.useRectifiedLinear = false; %instead of sigmoid
      this.params.Autoencoder.useRectifiedLog = false; %instead of sigmoid
      this.params.Autoencoder.useNoActFcn = false; %instead of sigmoid
      this.params.Autoencoder.useBinaryLinearInterp = false;
      this.params.Autoencoder.binaryLinearInterpScaling = 1;
      this.params.Autoencoder.binaryLinearInterpEMAconst = 0.99;
      
      this.params.Autoencoder.hiddenLinearBinaryThresholdUnitsFirst = false; %these are not considered when calculating gradients
      this.params.Autoencoder.hiddenSigmoidBinaryThresholdUnitsFirst = false; %these are not considered when calculating gradients

      this.params.Autoencoder.useSoftmax = false;
      this.params.Autoencoder.useNormMean = false;
      this.params.Autoencoder.normMeanEpsilon = 1e-6;
      
      this.params.Autoencoder.hiddenLinearBinaryThresholdUnits = false; %these are not considered when calculating gradients
      this.params.Autoencoder.hiddenSigmoidBinaryThresholdUnits = false; %these are not considered when calculating gradients
     
      this.params.Autoencoder.maskingNoiseFractionHidden = 0;        % fraction of inputs to set to 0
      this.params.Autoencoder.maskingNoiseFractionHiddenRescale = false;        % if inputs are rescaled to the variance before masking
      this.params.Autoencoder.maskingNoiseFractionHiddenEffectGradient = true; % if the masked y2 should be used when backpropagating errors
      
      this.params.Autoencoder.noBiasBack = false;

      this.params.Autoencoder.reconstrSigmoid = false;
      this.params.Autoencoder.reconstrRectifiedLinear = false;

      % OPTIMIZATION PARAMS:
      this.params.Autoencoder.alpha = 1;        % weight of Autoencoder 
      %WARNING: This is not scaled proportional to number units in a3. Therefore reduce alpha when increasing image size
      this.params.Autoencoder.useLinearAutoencError = false;
      this.params.Autoencoder.useCrossEntropyError = false;
      this.params.Autoencoder.beta = 5;         % weight of sparsity penalty term
      this.params.Autoencoder.sparsityParam = 0.035; % desired average activation of the hidden units.
      this.params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
      this.params.Autoencoder.gammaPvalue = 1;        % Psi_gamma = gamma * sqrt( NghMat * y2.^2 + epsilon ) ^ gammaPvalue
      this.params.Autoencoder.topoNghFcn = [];
      this.params.Autoencoder.topoNgh = [3 3];
      this.params.Autoencoder.topoGridDims = [20 20];
      this.params.Autoencoder.topoPeriodicBoundary = [true true];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      this.params.Autoencoder.mu = 0;        % weight of modular topographic sparsity penalty term
      this.params.Autoencoder.muPvalue = 1;        % Psi_mu = mu * ( muFeatureKernel * (muConvKernel * y2.^2) + muEpsilon ) ^ muPvalue
      this.params.Autoencoder.muConvKernel = ones(3,3);
      this.params.Autoencoder.muFeatureKernel = []; % if empty then no combination in feature dimension, otherwise it should be a matrix of size (hiddenSize,hiddenSize)
      this.params.Autoencoder.muEpsilon = 1e-2;
      this.params.Autoencoder.nu = 0;
      this.params.Autoencoder.nuEMAconst = 0.9;
      this.params.Autoencoder.nuCovWeightMat = [];
      this.params.Autoencoder.lambda = 3e-3;    % weight of L2-decay parameter
      this.params.Autoencoder.lambdaBackScale = 1;    % weight of L2-decay parameter scaling for backward conenctions
      this.params.Autoencoder.lambdaL1 = 0 ;    % weight of L1-decay parameter
      this.params.Autoencoder.lambdaL1BackScale = 1;    % weight of L1-decay parameter
      this.params.Autoencoder.fixL2WeightBack = false; %if true, remember to also set both backScale to 0
      this.params.Autoencoder.fixL2WeightBackPoolPerHiddenNeuron = false;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0 ; %if not zero then smooth it
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 0 ;
    
      % LOGS AND SAVES:
      this.params.Autoencoder.saveinterval = Inf;
      this.params.Autoencoder.saveintervalForwConn = Inf;
      this.params.Autoencoder.saveintervalBackConn = Inf;
      this.params.Autoencoder.saveintervalTrainInfo = Inf;
      this.params.Autoencoder.saveintervalLog = 10; %save to log variable
      this.params.Autoencoder.saveintervalpng = 10;
      this.params.Autoencoder.saveintervalpngKeep = 10;
      this.params.Autoencoder.savepng = false;
      this.params.Autoencoder.plotPermute = []; %function to convert W to 5 dim matrix.
      this.params.Autoencoder.reloadSaveinterval = 5;
      this.params.Autoencoder.debugOutput = false;
      
      % VALIDATION:
      this.params.Autoencoder.validationSetsize = 0;
      this.params.Autoencoder.validationInterval = Inf; %save to seperate file
      this.params.Autoencoder.validationIntervalLog = 10;
      this.params.Autoencoder.validationSetIds = [];
      this.params.Autoencoder.validationSetImageIds = []; % will only be used if also numberOfImagesPerPatchReload is set
      this.params.Autoencoder.validationSkipGrad = true;

      % ITERATION META-PARAMS:
      this.params.Autoencoder.resetRandstreamEachReload = false;
      this.params.Autoencoder.resetRandstreamEachIter = false;
      this.params.Autoencoder.resetRandstreamEachEval = false;
      this.params.Autoencoder.fixedBatches = false;
      this.params.Autoencoder.batchsize = []; %[] means full batch
      this.params.Autoencoder.numberOfPatchReloads = 1;
      this.params.Autoencoder.numberOfImagesPerPatchReload = [];
      this.params.Autoencoder.useMinFuncGrad = true;
      this.params.Autoencoder.initialItersWithoutLearning = 0;
      
      this.params.minFuncGrad.Method = 'rmsprop';
      this.params.minFuncGrad.maxIter = 500;
      this.params.minFuncGrad.learnrate = 1e-5;
      this.params.minFuncGrad.EMAconst = 0.1;
      this.params.minFuncGrad.momentum = 0.9; % = 1 - 2/(lambda*2.8854+1) for halflife lambda
      this.params.minFuncGrad.display = 'on';
      this.params.minFuncGrad.displayEvery = 1;
      
      this.params.minFunc.Method = 'lbfgs';
      this.params.minFunc.maxIter = 500;
      this.params.minFunc.display = 'on';
      this.params.minFunc.optTol = 1e-5; % set to 0 for minibatch
      this.params.minFunc.progTol = 1e-9; % set to -1 for minibatch
      this.params.minFunc.Damped = 0; %maybe set this to 1 for minibatch
      this.params.minFunc.Corr = 100; %maybe reduce this for minibatch ?
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      global logVars;
      global EMA_z2pos;
      global EMA_y2_pairs;
      global EMA_y2_square;
      global EMA_y2;
      
      if this.numJobs > 1
        this.tempsavepath = fullfile(this.temppath,num2str(this.currJobid));
        mkdir(this.tempsavepath);
      else
        this.tempsavepath = this.temppath;
      end
      
      inputfolder = fullfile(this.workpath, this.params.Autoencoder.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder, this.params.Autoencoder.inActFilenames );
      
      
      if ~isempty(this.params.Autoencoder.inputApplyConn)
        this.inputApplyConn = load(fullfile(this.workpath,this.params.Autoencoder.inputApplyConn));
      end
      
      
      if this.params.Autoencoder.validationSetsize~=0
        if ~isempty(this.params.Autoencoder.numberOfImagesPerPatchReload)
          if isempty(this.params.Autoencoder.validationSetImageIds)
            validationSetImageIds=randi(length(pathlist),1,this.params.Autoencoder.numberOfImagesPerPatchReload);
          else
            validationSetImageIds=this.params.Autoencoder.validationSetImageIds;
          end
          pathlistValid = pathlist(validationSetImageIds);
          filelistValid = filelist(validationSetImageIds);
          pathlist = pathlist(setdiff(1:length(pathlist),validationSetImageIds));
          filelist = filelist(setdiff(1:length(filelist),validationSetImageIds));
          validPatches = this.preparePatches(pathlistValid,filelistValid,this.params.Autoencoder.validationSetsize);
        else
          validPatches = this.preparePatches(pathlist,filelist,this.params.Autoencoder.validationSetsize);
        end
        
        validRandStream = RandStream('mt19937ar','Seed',randi(2^32-1));
        this.randstreamStateBeforeValidation = validRandStream.State;
      else
        validPatches = [];
      end
      
      this = this.setPatchsizes();
      
      if this.params.Autoencoder.continue && exist(fullfile(this.temppath,['currIterJob' num2str(this.currJobid) '.mat']),'file')
        tmp = load(fullfile(this.temppath,['currIterJob' num2str(this.currJobid) '.mat']));
        this.tempsavepath = tmp.this.tempsavepath;
        this.forwardRfSize = tmp.this.forwardRfSize;
        this.backwardRfSize = tmp.this.backwardRfSize;
        this.nghMat = tmp.this.nghMat;
        this.randstreamStateBeforeEval = tmp.this.randstreamStateBeforeEval;
        this.randstreamStateBeforeIter = tmp.this.randstreamStateBeforeIter;
        this.randstreamStateBeforeReload = tmp.this.randstreamStateBeforeReload;
        this.randstreamStateBeforeValidation = tmp.this.randstreamStateBeforeValidation;
        this.theta = tmp.this.theta;
        this.EMA_z2pos = tmp.this.EMA_z2pos;
        this.EMA_y2_pairs = tmp.this.EMA_y2_pairs;
        this.EMA_y2_square = tmp.this.EMA_y2_square;
        this.EMA_y2 = tmp.this.EMA_y2;
        this.nextReloadId = tmp.this.nextReloadId;
        this.log = tmp.this.log;
        this.gradStates = tmp.this.gradStates;
        clear tmp;
      else
        if this.params.Autoencoder.continueBatch
          inWpath = fullfile(this.workpath,this.params.Autoencoder.continueBatchInWeightsFolder);
          if this.numJobs > 1
            inWpath = fullfile(inWpath,num2str(this.currJobid));
          end
          tmp = load(fullfile(inWpath,this.params.Autoencoder.continueBatchInForwConnFilenames));
          % !!!!!!! not working correctly! TODO: reorder dimensions of W!
          W1 = reshape(tmp.W,this.dimW1);
          b1 = tmp.b;
          tmp = load(fullfile(inWpath,this.params.Autoencoder.continueBatchInBackConnFilenames));
          % !!!!!!! not working correctly! TODO: reorder dimensions of W!
          W2 = reshape(tmp.W,this.dimW2);
          b2 = tmp.b;
          clear tmp;
        else
          if this.params.Autoencoder.initWGaussian
            W1 = randn(this.dimW1) * this.params.Autoencoder.forwInitMaxWeight;
            W2 = randn(this.dimW2) * this.params.Autoencoder.backInitMaxWeight;
          else
            W1 = ( 2*rand(this.dimW1)-1 ) * this.params.Autoencoder.forwInitMaxWeight;
            W2 = ( 2*rand(this.dimW2)-1 ) * this.params.Autoencoder.backInitMaxWeight;
          end
          
          if this.params.Autoencoder.initWExponent ~= 1
            W1 = W1 .^ this.params.Autoencoder.initWExponent;
            W2 = W2 .^ this.params.Autoencoder.initWExponent;
          end

          if ~isempty(this.params.Autoencoder.forwInitScaleDistFcn)
            W1 = Autoencoder.scaleWeightWithDist(W1,this.params.Autoencoder.forwInitScaleDistFcn);
          end
          if ~isempty(this.params.Autoencoder.backInitScaleDistFcn)
            W2 = Autoencoder.scaleWeightWithDist(W2,this.params.Autoencoder.backInitScaleDistFcn);
          end

          b1 = zeros(this.dimW1(4:end));
          b2 = zeros([1 1 this.dimW2(6)]);

        end
        this.theta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
        
        if ~isempty(this.params.Autoencoder.loadPatchesAndTheta)
          load(this.params.Autoencoder.loadPatchesAndTheta)
        end

        if this.params.Autoencoder.gamma~=0
          if isempty(this.params.Autoencoder.topoNghFcn)
            this.nghMat = linearDecoderNghMat( this.params.Autoencoder.topoNgh,this.params.Autoencoder.topoGridDims,this.params.Autoencoder.topoPeriodicBoundary);
          else
            this.nghMat = linearDecoderNghWinMat( this.params.Autoencoder.topoNgh,this.params.Autoencoder.topoGridDims,this.params.Autoencoder.topoPeriodicBoundary, this.params.Autoencoder.topoNghFcn);
          end
        end
        
        this.selectMinibatch();
      
        globStream = Autoencoder.getGlobRng();
        if this.params.Autoencoder.resetRandstreamEachReload
          this.randstreamStateBeforeReload = globStream.State;
        end

        if ~isfield(this.log,'train')
          this.log.train = [];
          
        end
        if ~isfield(this.log,'validation')
          this.log.validation = [];
        end
        
      end
      
      if isempty(EMA_z2pos)
        EMA_z2pos = 0.01*ones(this.dimW1(4:6));
      end
      
      
      logVars = this.log;
      
      for reloadId = this.nextReloadId:this.params.Autoencoder.numberOfPatchReloads-1
        
        numTrainPatches = this.params.Autoencoder.inSamplesDims(4) - this.params.Autoencoder.validationSetsize;
        if ~isempty(this.params.Autoencoder.numberOfImagesPerPatchReload)
          imgFileIds = randi(length(pathlist),1,this.params.Autoencoder.numberOfImagesPerPatchReload);
          tic;
          patches = this.preparePatches(pathlist(imgFileIds),filelist(imgFileIds),numTrainPatches);
          if this.params.Autoencoder.debugOutput
            disp(['Loading time for new images: ' num2str(toc)]);
          end
        else
          patches = this.preparePatches(pathlist,filelist,numTrainPatches);
        end
        
        if this.params.Autoencoder.resetRandstreamEachReload
          globStream.State = this.randstreamStateBeforeReload;
        end
        if this.params.Autoencoder.resetRandstreamEachIter
          this.randstreamStateBeforeIter = globStream.State;
        end
        if this.params.Autoencoder.resetRandstreamEachEval
          this.randstreamStateBeforeEval = globStream.State;
        end
        
        if length(this.params.Autoencoder.hiddenSize)>1
          old = this.currentNumberOfUnits;
          this.currentNumberOfUnits = this.params.Autoencoder.hiddenSize(min(length(this.params.Autoencoder.hiddenSize),reloadId+1));
          if old~=this.currentNumberOfUnits
            %reset grad states because length is changing
            this.gradStates = [];
            this.tiledConvF = TiledConv(this.tiledConvF.tileSizeX,this.tiledConvF.tileSizeY,this.tiledConvF.fIn,this.currentNumberOfUnits);
            this.tiledConvB = TiledConv(this.tiledConvB.tileSizeX,this.tiledConvB.tileSizeY,this.currentNumberOfUnits,this.tiledConvF.fOut);
            this.tiledConvInv = TiledConv(this.tiledConvB.tileSizeX,this.tiledConvB.tileSizeY,this.tiledConvF.fOut,this.currentNumberOfUnits);
          end
          [W1, W2, b1, b2] = this.decodeTheta( this.theta );
          W1 = W1(:,:,:,:,:,1:this.currentNumberOfUnits);
          W2 = W2(:,:,1:this.currentNumberOfUnits,:,:,:);
          b1 = b1(:,:,1:this.currentNumberOfUnits);
          tempTheta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
        else
          this.currentNumberOfUnits = this.dimW1(6);%this.params.Autoencoder.hiddenSize;
          tempTheta = this.theta;
        end
        
        if reloadId==this.nextReloadId && this.params.Autoencoder.initialItersWithoutLearning
          for initIter=1:this.params.Autoencoder.initialItersWithoutLearning
            this.autoencoderCost(tempTheta, patches, false, true, false);
          end
        end
        
        if this.params.Autoencoder.useMinFuncGrad
          options = this.params.minFuncGrad;
          options.outputFcn = @(optTheta,iterationType,i,funEvals,f,t,gtd,g,d,optCond) this.iterFcnSave( optTheta,iterationType,i+reloadId*options.maxIter,funEvals,f,t,gtd,g,d,optCond, validPatches);
          %         save(fullfile(this.tempsavepath,'initialParams.mat'),'-v7.3','theta','this','patches','options');%
          [tempTheta, ~,~,~,this.gradStates] = minFuncGrad( @(p) autoencoderCost(this, p, patches, false, true, false), tempTheta, options, this.gradStates);
        else
          options = this.params.minFunc;
          options.outputFcn = @(optTheta,iterationType,i,funEvals,f,t,gtd,g,d,optCond) this.iterFcnSave( optTheta,iterationType,i+reloadId*options.maxIter,funEvals,f,t,gtd,g,d,optCond, validPatches);
          %         save(fullfile(this.tempsavepath,'initialParams.mat'),'-v7.3','theta','this','patches','options');
          tempTheta = minFunc( @(p) autoencoderCost(this, p, patches, false, true, false), tempTheta, options);
        end
        
        if length(this.params.Autoencoder.hiddenSize)>1
          % decode theta and return to full size weight matrices:
          [W1, W2, b1, b2] = this.decodeTheta( tempTheta );
          
          % reencode the full weight matrices:
          this.theta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
        else
          this.theta = tempTheta;
        end
        
        % save this class object:
        if mod(reloadId,this.params.Autoencoder.reloadSaveinterval)==0
          this.nextReloadId = reloadId + 1;
          this.log = logVars;
          if this.params.Autoencoder.useBinaryLinearInterp
            this.EMA_z2pos = EMA_z2pos;
          end
          if this.params.Autoencoder.nu
            this.EMA_y2_pairs = EMA_y2_pairs;
            this.EMA_y2_square = EMA_y2_square;
            this.EMA_y2 = EMA_y2;
          end
          save(fullfile(this.temppath,['currIterJob' num2str(this.currJobid) '.mat']),'this');
        end
      end
      
      %% save results:
      savepath = fullfile(this.workpath,this.params.Autoencoder.outWeightsFolder);
      if this.numJobs > 1
        savepath = fullfile(savepath,num2str(this.currJobid));
      end
      mkdir(savepath);
      
      [forwConn, backConn] = this.createConn(this.theta); %#ok<ASGLU,NASGU>
      save(fullfile(savepath,'forwConn.mat'),'-struct','forwConn')
      save(fullfile(savepath,'backConn.mat'),'-struct','backConn')
      
      disp('Finished saving connections');
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function [this] = setPatchsizes(this)
      global EMA_y2_pairs;
      global EMA_y2_square;
      global EMA_y2;
      
      patchF = this.params.Autoencoder.patchDimForward;
      patchB = this.params.Autoencoder.patchDimBackward;
      fIn = this.params.Autoencoder.inSamplesDims(3)*this.params.Autoencoder.inputSubsampling.^2;
      fOut = this.params.Autoencoder.hiddenSize(end);
      if this.params.Autoencoder.useTiledConv
        tileSizeF=patchF-1;
        tileSizeB=patchB-1;
        this.tiledConvF = TiledConv(tileSizeF,tileSizeF,fIn,fOut);
        this.tiledConvB = TiledConv(tileSizeB,tileSizeB,fOut,fIn);
        this.tiledConvBInv = TiledConv(tileSizeB,tileSizeB,fIn,fOut);
      else
        tileSizeF=1;
        tileSizeB=1;
      end
      this.forwardRfSize = [patchF patchF fIn];
      this.backwardRfSize = [patchB patchB fOut];
      this.dimW1 = [patchF patchF fIn tileSizeF tileSizeF fOut];
      this.dimW2 = [patchB patchB fOut tileSizeB tileSizeB fIn];
      
      if this.params.Autoencoder.nu
        EMA_y2_pairs = zeros(prod(this.dimW1(4:6)),prod(this.dimW1(4:6)));
        EMA_y2_square = zeros(1,prod(this.dimW1(4:6)));
        EMA_y2 = zeros(1,prod(this.dimW1(4:6)));
      end
      
    end
    
    function [patches] = preparePatches(this,pathlist,filelist,numSamples)
      
      if this.params.Autoencoder.inputApplyConn
        patches = loadRandomPatches( ...
          this.params.Autoencoder.inSamplesDims(1), ...
          this.params.Autoencoder.inSamplesDims(2),...
          this.params.Autoencoder.inSamplesDims(3), ...
          numSamples, ...
          pathlist, ...
          filelist, ...
          this.params.Autoencoder.inSamplesBorderBuffer,...
          this.inputApplyConn);
      else
        patches = loadRandomPatches( ...
          this.params.Autoencoder.inSamplesDims(1), ...
          this.params.Autoencoder.inSamplesDims(2),...
          this.params.Autoencoder.inSamplesDims(3), ...
          numSamples, ...
          pathlist, ...
          filelist, ...
          this.params.Autoencoder.inSamplesBorderBuffer );
      end
      
      if this.params.Autoencoder.inputSubsampling > 1
        % now collapse the subsampled dimensions into the feature dim:
        patches = reshape(patches, [...
          this.params.Autoencoder.inputSubsampling ...
          size(patches,1)/this.params.Autoencoder.inputSubsampling ...
          this.params.Autoencoder.inputSubsampling ...
          size(patches,2)/this.params.Autoencoder.inputSubsampling ...
          size(patches,3) ...
          size(patches,4)]);
        patches = permute(patches,[2 4 1 3 5 6]);
        patches = reshape(patches,[...
          size(patches,1) ...
          size(patches,2) ...
          size(patches,3)*size(patches,4)*size(patches,5) ...
          size(patches,6)]);
      end
    end
    
    function [forwConn, backConn] = createConn(this, theta)
      global EMA_z2pos;
      
      [forwConn.W, backConn.W, forwConn.b, backConn.b] = this.decodeTheta( theta );
      %W has dims [dxin dyin fin dxout dyout fout]
      %b has dims [dxout dyout fout]
      
      if this.params.Autoencoder.inputSubsampling > 1
        
        % inputSubsampling is so far all collapsed in fin
        % so here it is extracted and collapsed into dxin and dyin
        forwConn.W = reshape(forwConn.W,[...
          size(forwConn.W,1)...
          size(forwConn.W,2)...
          this.params.Autoencoder.inputSubsampling...
          this.params.Autoencoder.inputSubsampling...
          size(forwConn.W,3)/this.params.Autoencoder.inputSubsampling^2 ...
          size(forwConn.W,4)...
          size(forwConn.W,5)...
          size(forwConn.W,6)]);
        forwConn.W = permute(forwConn.W,[3 1 4 2 5 6 7 8]);
        forwConn.W = reshape(forwConn.W,[...
          size(forwConn.W,1)*size(forwConn.W,2)...
          size(forwConn.W,3)*size(forwConn.W,4)...
          size(forwConn.W,5)...
          size(forwConn.W,6)...
          size(forwConn.W,7)...
          size(forwConn.W,8)]);
        
        %now reshape backward weights:
        backConn.W = reshape(backConn.W,[...
          size(backConn.W,1)... %dxin
          size(backConn.W,2)... %dyin
          size(backConn.W,3)... %fin
          size(backConn.W,4)... %tilexout
          size(backConn.W,5)... %tileyout
          this.params.Autoencoder.inputSubsampling... %dxout subsampling
          this.params.Autoencoder.inputSubsampling... %dyout subsampling
          size(backConn.W,6)/this.params.Autoencoder.inputSubsampling^2]); %fout
        backConn.W = permute(backConn.W,[1 2 3 6 4 7 5 8]);
        backConn.W = reshape(backConn.W,[...
          size(backConn.W,1)... %dxin
          size(backConn.W,2)... %dyin
          size(backConn.W,3)... %fin
          size(backConn.W,4)*size(backConn.W,5)... %dxout*tilexout
          size(backConn.W,6)*size(backConn.W,7)... %dyout*tileyout
          size(backConn.W,8)]); %fout
      end
      
      % save forward connections:
      forwConn.inputSubsampling = this.params.Autoencoder.inputSubsampling;
      forwConn.shiftOutputdims = false; %at the moment no tiled conv implemented
      
      if this.params.Autoencoder.useRectifiedLinear
        forwConn.actFcn = @(x) max(x,0);
      elseif this.params.Autoencoder.useRectifiedLog
        forwConn.actFcn = @(x) log(max(x,0)+1); %#ok<*CPROP>
      elseif this.params.Autoencoder.useNoActFcn
        forwConn.actFcn = [];
      elseif this.params.Autoencoder.useBinaryLinearInterp
        if ~isempty(EMA_z2pos)
          forwConn.W = bsxfun(@rdivide,forwConn.W, this.params.Autoencoder.binaryLinearInterpScaling*shiftdim(EMA_z2pos,-3));
          forwConn.b = bsxfun(@rdivide,forwConn.b, this.params.Autoencoder.binaryLinearInterpScaling*EMA_z2pos);
        end
        forwConn.actFcn = @(x) min(max(x,0),1);
      else
        forwConn.actFcn = @(x) 1 ./ (1 + exp(-x));
      end
      
      if this.params.Autoencoder.useSoftmax
        forwConn.actFcn2 = @(x) feval(@(expx) bsxfun(@rdivide,expx,sum(expx,3)), exp(x));
      elseif this.params.Autoencoder.useNormMean
        forwConn.actFcn2 = @(x) bsxfun(@rdivide,x,sum(x,3)+this.params.Autoencoder.normMeanEpsilon);
      end
      
      % save backward connections:
      backConn.shiftOutputdims = false;
      backConn.actFcn = []; % at the moment only linear
      
    end
    
  end
  
  methods (Static)
    
    function W = scaleWeightWithDist(W,distInitFcn)
      distMat = sqrt(bsxfun(@plus,((1:size(W,1))-(size(W,1)+1)/2).^2', ((1:size(W,2))-(size(W,2)+1)/2).^2));
      distMat = feval(distInitFcn,distMat);
      distMat = distMat / mean(distMat(:));
      W=bsxfun(@times,W,distMat);
    end
  end
  
  methods
    
    function selectMinibatch( this )
      global minibatchIds;
      if ~isempty(this.params.Autoencoder.batchsize)
        if this.params.Autoencoder.fixedBatches
          startId=randi(this.params.Autoencoder.inSamplesDims(4)-this.params.Autoencoder.validationSetsize-this.params.Autoencoder.batchsize+1,1,1);
          minibatchIds=startId:(startId+this.params.Autoencoder.batchsize-1);
        else
          tmp = randperm(this.params.Autoencoder.inSamplesDims(4)-this.params.Autoencoder.validationSetsize);
          minibatchIds=tmp(1:this.params.Autoencoder.batchsize);
        end
      end
    end
    
    
    function [stop] = iterFcnSave( this, optTheta, iterationType,i,funEvals,f,t,gtd,g,d,optCond, validPatches ) %#ok<INUSL>
      global minibatchIds;
      global logVars;
      
      if mod(i,this.params.Autoencoder.saveintervalpng)==0 || mod(i,this.params.Autoencoder.saveintervalpngKeep)==0 || mod(i,this.params.Autoencoder.saveinterval)==0
        [forwConn, backConn] = this.createConn(optTheta); %#ok<NASGU>
      end
      
      if mod(i,this.params.Autoencoder.saveintervalForwConn)==0
        save(fullfile(this.tempsavepath,['forwConnIter' num2str(i) '.mat']),'-struct','forwConn');
      end
      if mod(i,this.params.Autoencoder.saveintervalBackConn)==0
        save(fullfile(this.tempsavepath,['backConnIter' num2str(i) '.mat']),'-struct','backConn');
      end
      if mod(i,this.params.Autoencoder.saveintervalTrainInfo)==0
        save(fullfile(this.tempsavepath,['trainInfoIter' num2str(i) '.mat']),'iterationType','i','funEvals','f','t','gtd','d','optCond');
      end
      
      if mod(i,this.params.Autoencoder.saveinterval)==0
        save(fullfile(this.tempsavepath,['forwConnIter' num2str(i) '.mat']),'-struct','forwConn');
        save(fullfile(this.tempsavepath,['backConnIter' num2str(i) '.mat']),'-struct','backConn');
        save(fullfile(this.tempsavepath,['trainInfoIter' num2str(i) '.mat']),'iterationType','i','funEvals','f','t','gtd','d','optCond');
      end
      
      if mod(i,this.params.Autoencoder.saveintervalLog)==0
        logVars.train(end+1).iterationType = iterationType;
        logVars.train(end).i = i;
        logVars.train(end).funEvals = funEvals;
        logVars.train(end).f = f;
        logVars.train(end).optCond = optCond;
      end
      
      if mod(i,this.params.Autoencoder.saveintervalpng)==0 || mod(i,this.params.Autoencoder.saveintervalpngKeep)==0
        if this.params.Autoencoder.savepng
          pngpathKeep = fullfile(this.tempsavepath,['rfsForwConnIter' num2str(i) '.png']);
          pngpath = fullfile(this.temppath,['rfsForwConnLastIterJob' num2str(this.currJobid) '.png']);
          this.plotRfs( pngpath, forwConn );
          if mod(i,this.params.Autoencoder.saveintervalpngKeep)==0
            copyfile(pngpath,pngpathKeep,'f');
          end
        end
      end
      
      
      %% do validation set:
      if this.params.Autoencoder.validationSetsize~=0
        if mod(i,this.params.Autoencoder.validationInterval)==0 || mod(i,this.params.Autoencoder.validationIntervalLog)==0
          
          % set both randstream states such that it will be used during validation:
          globStream = Autoencoder.getGlobRng();
          tmpGlobalRandstreamState = globStream.State;
          globStream.State = this.randstreamStateBeforeValidation;
          this.randstreamStateBeforeEval = this.randstreamStateBeforeValidation;
          
          % change minibatch for validation
          if ~isempty(this.params.Autoencoder.batchsize)
            tmp = minibatchIds;
            minibatchIds = 1:this.params.Autoencoder.validationSetsize;
          end
          
          % run validation set:
          iterationType = struct();
          tic;
          
          [cost,~,optout] = this.autoencoderCost( optTheta, validPatches, this.params.Autoencoder.validationSkipGrad, false, true);
            
          iterationType.optout = optout;
          disp(['iter ' num2str(i) ' validation cost: ' num2str(cost) ' validation time: ' num2str(toc)])
          iterationType.cost = cost;
          iterationType.iteration = i;
          
          % save validation into log variable
          if mod(i,this.params.Autoencoder.validationIntervalLog)==0
            logVars.validation(end+1).iterationType = iterationType;
            logVars.validation(end).i = i;
            logVars.validation(end).cost = cost;
          end
          
          % save validation to seperate file
          if mod(i,this.params.Autoencoder.validationInterval)==0
            save(fullfile(this.tempsavepath,['validInfoIter' num2str(i) '.mat']),'iterationType','i','funEvals','f','t','gtd','d','optCond');          
          end
          
          % reset minibatch to previous values
          if ~isempty(this.params.Autoencoder.batchsize)
            minibatchIds = tmp;
          end
          
          % reset randstream state back to the state before validation:
          globStream.State = tmpGlobalRandstreamState;
          
        end
      end
      
      this.selectMinibatch();
      
      if this.params.Autoencoder.resetRandstreamEachIter
        globStream = Autoencoder.getGlobRng();
        globStream.State = this.randstreamStateBeforeIter;
      end
      
      stop=false;
    end
    
    function plotRfs( this, pngpath, forwConn )
      if nargin<3
        forwConn = this.createConn(this.theta);
      end
      
      totalNumFeaturesForw = size(forwConn.W,4)*size(forwConn.W,5)*size(forwConn.W,6);
      plotSizeX = round(sqrt(totalNumFeaturesForw));
      while mod(totalNumFeaturesForw/plotSizeX,1)~=0 && plotSizeX>1
        plotSizeX = plotSizeX - 1;
      end
      
      if this.params.Autoencoder.inputApplyConn
        
        dimPrevLayer = size(this.inputApplyConn.W);
        dimThisLayer = size(forwConn.W);
        
        convRfs = zeros(dimPrevLayer(1)+dimThisLayer(1)-1,dimPrevLayer(2)+dimThisLayer(2)-1,dimPrevLayer(3),dimThisLayer(6));
        for i=1:dimPrevLayer(6)
          for j=1:dimThisLayer(6)
            for k=1:dimPrevLayer(3)
              convRfs(:,:,k,j) = convRfs(:,:,k,j) + conv2(forwConn.W(:,:,i,1,1,j), this.inputApplyConn.W(:,:,k,1,1,i));
            end
          end
        end
        plotColorFeatures( reshape(convRfs,[size(convRfs,1) size(convRfs,2) size(convRfs,3) plotSizeX totalNumFeaturesForw/plotSizeX]), true, pngpath, true );
      else
        if size(forwConn.W,3)==6
          % if we have a concatinated input of two times RGB features (one for means and one for edges)
          plotColorFeatures( reshape(forwConn.W(:,:,1:3,:,:,:),[size(forwConn.W,1) size(forwConn.W,2) 3 plotSizeX totalNumFeaturesForw/plotSizeX]), false, [pngpath 'Cat1'], true );
          plotColorFeatures( reshape(forwConn.W(:,:,4:6,:,:,:),[size(forwConn.W,1) size(forwConn.W,2) 3 plotSizeX totalNumFeaturesForw/plotSizeX]), false, [pngpath 'Cat2'], true );
        else
          if ~isempty(this.params.Autoencoder.plotPermute)
            plotColorFeatures( feval(this.params.Autoencoder.plotPermute,forwConn.W), true, pngpath, true );
          else
            if this.params.Autoencoder.useTiledConv
              tileSize = totalNumFeaturesForw / this.params.Autoencoder.hiddenSize;
              plotColorFeatures( reshape(forwConn.W,[size(forwConn.W,1) size(forwConn.W,2) size(forwConn.W,3) sqrt(tileSize) sqrt(tileSize)*this.params.Autoencoder.hiddenSize]), true, pngpath, true );
            else
              plotColorFeatures( reshape(forwConn.W,[size(forwConn.W,1) size(forwConn.W,2) size(forwConn.W,3) plotSizeX totalNumFeaturesForw/plotSizeX]), true, pngpath, true );
            end
          end
        end
      end
      
    end
    
    function [W1, W2, b1, b2] = decodeTheta( this, theta, currentNumberOfUnits )
      
      if nargin<3
%         currentNumberOfUnits = this.backwardRfSize(3);
        currentNumberOfUnits = this.dimW1(6);
        
      end
      
      dimW1 = this.dimW1;
      dimW2 = this.dimW2;
      
      dimW1(6) = currentNumberOfUnits;
      
%       dimW1 = [this.forwardRfSize(1), this.forwardRfSize(2), this.forwardRfSize(3), 1, 1, currentNumberOfUnits];
%       dimW2 = [this.backwardRfSize(1), this.backwardRfSize(2), currentNumberOfUnits, 1, 1, this.forwardRfSize(3)];
      dimb1 = dimW1(4:6);%[1, 1, currentNumberOfUnits];
      dimb2 = [1 1 dimW2(6)];%[1, 1, this.forwardRfSize(3)];
      
      if prod(dimW1)+prod(dimW2)+prod(dimb1)+prod(dimb2) ~= numel(theta)
        
        %full weight matrices:
        [W1, W2, b1, b2] = this.decodeTheta( this.theta );
        
        %current part weight matrices:
        [W1cur, W2cur, b1cur, b2cur] = this.decodeTheta( theta, this.currentNumberOfUnits );
        
        W1(1:size(W1cur,1),1:size(W1cur,2),1:size(W1cur,3),1:size(W1cur,4),1:size(W1cur,5),1:size(W1cur,6)) = W1cur;
        W2(1:size(W2cur,1),1:size(W2cur,2),1:size(W2cur,3),1:size(W2cur,4),1:size(W2cur,5),1:size(W2cur,6)) = W2cur;
        b1(1:size(b1cur,1),1:size(b1cur,2),1:size(b1cur,3)) = b1cur;
        b2(1:size(b2cur,1),1:size(b2cur,2),1:size(b2cur,3)) = b2cur;
        
      else
      
        W1 = theta(1:prod(dimW1));
        W2 = theta(prod(dimW1)+1:prod(dimW1)+prod(dimW2));
        b1 = theta(prod(dimW1)+prod(dimW2)+1:prod(dimW1)+prod(dimW2)+prod(dimb1));
        b2 = theta(prod(dimW1)+prod(dimW2)+prod(dimb1)+1:prod(dimW1)+prod(dimW2)+prod(dimb1)+prod(dimb2));

        W1 = reshape(W1,dimW1);
        W2 = reshape(W2,dimW2);
        b1 = reshape(b1,dimb1);
        b2 = reshape(b2,dimb2);
      
      end
      
    end
    
    
    %% This is the cost function:
    function [cost,grad,optout] = autoencoderCost(this, theta, data, skipGrad, noOptout, noStateChange)
      global minibatchIds;
      global EMA_z2pos;
      global EMA_y2_pairs;
      global EMA_y2_square;
      global EMA_y2;
      
      if nargin<4
        skipGrad = false;
      end
      if nargin<5 || isempty(noOptout)
        if nargout<3
          noOptout = true;
        else
          noOptout = false;
        end
      end
      if nargin<6
        noStateChange = false;
      end
      optout = struct();
      
      if noStateChange
        temp_minibatchIds = minibatchIds;
        temp_EMA_z2pos =  EMA_z2pos;
        temp_EMA_y2_pairs =  EMA_y2_pairs;
        temp_EMA_y2_square =  EMA_y2_square;
        temp_EMA_y2 =  EMA_y2;
      end
      
      if this.params.Autoencoder.resetRandstreamEachEval
        globStream = Autoencoder.getGlobRng();
        globStream.State = this.randstreamStateBeforeEval;
      end
      
%       if this.params.Autoencoder.debugOutput
%         tmpglobStream = Autoencoder.getGlobRng();
%         tmpglobStreamInitState = tmpglobStream.State;
%         disp(['globstream state: ' num2str(tmpglobStreamInitState(1:4)')])
%       end
      
      if ~isempty(this.params.Autoencoder.batchsize)
        data=data(:,:,:,minibatchIds);
      end
      params = this.params.Autoencoder;
      
      if ~isempty(this.currentNumberOfUnits)
        [W1, W2, b1, b2] = this.decodeTheta( theta, this.currentNumberOfUnits );
      else
        [W1, W2, b1, b2] = this.decodeTheta( theta );
      end
      % Now: W has dims [dxin dyin fin dxout dyout fout]
      % Now: b has dims [dxout dyout fout]
      
      W1 = reshape(W1,[size(W1,1) size(W1,2) size(W1,3) size(W1,4)*size(W1,5)*size(W1,6)]);
      W2 = reshape(W2,[size(W2,1) size(W2,2) size(W2,3) size(W2,4)*size(W2,5)*size(W2,6)]);
      b1 = reshape(b1,[size(b1,1)*size(b1,2)*size(b1,3) 1]);
      b2 = reshape(b2,[size(b2,1)*size(b2,2)*size(b2,3) 1]);
      
%       W1 = permute(W1,[1 2 3 6 4 5]);
%       W2 = permute(W2,[1 2 3 6 4 5]);
%       b1 = permute(b1,[3 1 2]);
%       b2 = permute(b2,[3 1 2]);
      
      % Now: W has dims [dxin dyin fin fout 1 1]
      % Now: b has dims [fout 1 1]
      % Remember that inputSubsampling is all collapsed in fin

      
      
      if params.fixL2WeightBack
        W2tilde = W2;
        if params.fixL2WeightBackPoolPerHiddenNeuron
          W2tildeNormalizer = sqrt( sum(sum(sum( W2tilde.^2, 1 ), 2) , 4) );
        else
          W2tildeNormalizer = sqrt( sum(sum(sum( W2tilde.^2, 1 ), 2) , 3) );
        end
        W2 = bsxfun(@rdivide, W2tilde, W2tildeNormalizer );
      end
        
      m=size(data,4);
      
      if params.beta ~= 0
        
        if params.maskingNoiseFraction || params.gaussianNoiseSigmaInput || params.maskingNoiseFractionHidden
          globStream = Autoencoder.getGlobRng();
          globStreamInitState = globStream.State;
        end
        
        avgAct2counter = 0;
        avgAct2=zeros([1 1 numel(b1)]);
        for k=1:m
          curData = data(:,:,:,k);
          [~,~,~,y2] = this.forwardPass(curData,W1,b1);
          avgAct2 = avgAct2 + sum(sum(y2,2),1);
          avgAct2counter = avgAct2counter + size(y2,1)*size(y2,2);
        end
        avgAct2 = avgAct2/avgAct2counter;
        
        if params.maskingNoiseFraction || params.gaussianNoiseSigmaInput || params.maskingNoiseFractionHidden
          globStream.State = globStreamInitState;
        end
      end
      
      costAutoenc = 0;
      costRegularize = (params.lambda / 2) * (sum(W1(:).^2) + params.lambdaBackScale*sum(W2(:).^2));
      if params.scaleLambdaL1WithDistExponent
        distMatW1 = sqrt(bsxfun(@plus,((1:size(W1,1))-(size(W1,1)+1)/2).^2', ((1:size(W1,2))-(size(W1,2)+1)/2).^2));
        distMatW1 = distMatW1.^params.scaleLambdaL1WithDistExponent;
        distMatW1 = distMatW1 / mean(distMatW1(:));
        
        distMatW2 = sqrt(bsxfun(@plus,((1:size(W2,1))-(size(W2,1)+1)/2).^2', ((1:size(W2,2))-(size(W2,2)+1)/2).^2));
        distMatW2 = distMatW2.^params.scaleLambdaL1WithDistExponent;
        distMatW2 = distMatW2 / mean(distMatW2(:));
      end
      
      if params.lambdaL1
        if params.smoothLambdaL1WithEpsilon>0
          W1L1 = sqrt(params.smoothLambdaL1WithEpsilon + W1.^2);
          W2L1 = sqrt(params.smoothLambdaL1WithEpsilon + W2.^2);
        else
          W1L1 = abs(W1);
          W2L1 = abs(W2);
        end
        
%           if params.scaleLambdaL1WithDistExponent
%             costRegularizeL1 = 0;
%             costRegularizeL1 = costRegularizeL1 + params.lambdaL1 * sum(sqrt(params.smoothLambdaL1WithEpsilon + (distMatW1(:) .* sum(reshape(abs(W1),[size(W1,1)*size(W1,2) size(W1,3)*size(W1,4)*size(W1,5)*size(W1,5)]),2)).^2));
%             costRegularizeL1 = costRegularizeL1 + params.lambdaL1 * sum(sqrt(params.smoothLambdaL1WithEpsilon + (distMatW2(:) .* sum(reshape(abs(W2),[size(W2,1)*size(W2,2) size(W2,3)*size(W2,4)*size(W2,5)*size(W2,6)]),2)).^2));
%           else
%             costRegularizeL1 = params.lambdaL1 * (sum(sqrt(params.smoothLambdaL1WithEpsilon + W1(:).^2)) + sum(sqrt(params.smoothLambdaL1WithEpsilon + W2(:).^2)));
%           end
%         else
          if params.scaleLambdaL1WithDistExponent
            costRegularizeL1 = 0;
            costRegularizeL1 = costRegularizeL1 + params.lambdaL1 * sum(distMatW1(:) .* sum(reshape(W1L1,[size(W1,1)*size(W1,2) size(W1,3)*size(W1,4)*size(W1,5)*size(W1,6)]),2));
            costRegularizeL1 = costRegularizeL1 + (params.lambdaL1 * params.lambdaL1BackScale) * sum(distMatW2(:) .* sum(reshape(W2L1,[size(W2,1)*size(W2,2) size(W2,3)*size(W2,4)*size(W2,5)*size(W2,6)]),2));
          else
            costRegularizeL1 = params.lambdaL1 * (sum(W1L1(:)) + params.lambdaL1BackScale * sum(W2L1(:)));
          end
%         end
      else
        costRegularizeL1 = 0;
      end
      
      if params.beta~=0
        costSparseness = params.beta*sum( params.sparsityParam*log(params.sparsityParam./avgAct2(:)) + (1-params.sparsityParam)*log((1-params.sparsityParam)./(1-avgAct2(:))) );
      else
        costSparseness=0;
      end
      
      W1grad = zeros(size(W1));
      W2grad = zeros(size(W2));
      b1grad = zeros(size(b1));
      b2grad = zeros(size(b2));
      
      numUnits2 = (size(data,1)-this.forwardRfSize(1)+1)*(size(data,2)-this.forwardRfSize(2)+1);
      
      avgPsi_sparse=zeros([1 numel(b1)]);
      cost_mu = 0;
      cost_nu = 0;
      
      if ~noOptout
        sumSquaredDelta2Autoenc = 0;
        sumSquaredDelta2Sparse = 0;
        sumSquaredDelta2Toposparse = 0;
        sumSquaredDelta2mu = 0;
        sumDelta2Autoenc = 0;
        sumDelta2Sparse = 0;
        sumDelta2Toposparse = 0;
        sumDelta2mu = 0;
        
        sumActivation_z2 = zeros([1 1 numel(b1)]);
        sumActivationSquared_z2 = zeros([1 1 numel(b1)]);
        sumActivationFourthMoment_z2 = zeros([1 1 numel(b1)]);
        numActive_z2 = zeros([1 1 numel(b1)]);
        sumActivation_a2 = zeros([1 1 numel(b1)]);
        sumActivationSquared_a2 = zeros([1 1 numel(b1)]);
        sumActivationFourthMoment_a2 = zeros([1 1 numel(b1)]);
        numActive_a2 = zeros([1 1 numel(b1)]);
        sumActivation_y2 = zeros([1 1 numel(b1)]);
        sumActivationSquared_y2 = zeros([1 1 numel(b1)]);
        sumActivationFourthMoment_y2 = zeros([1 1 numel(b1)]);
        numActive_y2 = zeros([1 1 numel(b1)]);
        numPerFeatureType_layer2 = zeros([1 1 numel(b1)]);

      end
      
      if params.nu~=0
        d_EMA_y2_pairs = zeros(prod(this.dimW1(4:6)),prod(this.dimW1(4:6)));
        d_EMA_y2 = zeros(1,prod(this.dimW1(4:6)));
        d_EMA_y2_square = zeros(1,prod(this.dimW1(4:6)));
      end
      
      
      if params.nu~=0
        [cov_y2, var_y2, var_pair] = this.getVariances();
        
        delta_cov_y2 = params.nu * 2 * cov_y2 ./ var_pair;
        if ~isempty(params.nuCovWeightMat)
          delta_cov_y2 = delta_cov_y2 .* params.nuCovWeightMat;
        end
        delta_cov_y2(logical(eye(size(delta_cov_y2)))) = 0; %set diagonal to zero
        if ~isempty(params.nuCovWeightMat)
          delta_var_y2 = - params.nu * 2 * sum( bsxfun(@rdivide, params.nuCovWeightMat .* cov_y2.^2, var_y2), 2 )' ./ var_y2.^2;
        else
          delta_var_y2 = - params.nu * 2 * sum( bsxfun(@rdivide, cov_y2.^2, var_y2), 2 )' ./ var_y2.^2;
        end
        
        delta_EMA_y2_square = delta_var_y2;
        delta_EMA_y2 = 2*sum( bsxfun(@times, delta_cov_y2, -EMA_y2), 2)' - 2*delta_var_y2.*EMA_y2;%TODO: sum over 1 or 2nd dim?
        delta_EMA_y2_pairs = delta_cov_y2;
      end
      
      for k=1:m
        curData = data(:,:,:,k);
        
        [curData,z2,a2,y2,localNormalizer] = this.forwardPass(curData,W1,b1);
        if params.useTiledConv
          y2NumTilesX = size(curData,1)/this.tiledConvF.tileSizeX - 1;
          y2NumTilesY = size(curData,2)/this.tiledConvF.tileSizeY - 1;
        end
        
        if params.useBinaryLinearInterp
          mean_z2pos = sum(sum(max(z2,0),1),2) ./ sum(sum(z2>0, 1),2);
          mean_z2pos(isnan(mean_z2pos)) = EMA_z2pos(isnan(mean_z2pos));
          EMA_z2pos = EMA_z2pos * params.binaryLinearInterpEMAconst + reshape(mean_z2pos,size(EMA_z2pos)) * (1-params.binaryLinearInterpEMAconst);
        end
        
        if params.gamma~=0
          y2reduced = reshape(y2,[size(y2,1)*size(y2,2) size(y2,3)]);
          Psi_sparse = (y2reduced.^2 * this.nghMat + params.topoEpsilon).^(0.5 * params.gammaPvalue);
          avgPsi_sparse = avgPsi_sparse + sum(Psi_sparse,1);
        end
        
        if params.mu~=0
          y2SqrTiled = reshape(y2.^2,[y2NumTilesX y2NumTilesY this.tiledConvF.tileSizeX this.tiledConvF.tileSizeY this.tiledConvF.fOut]);
          y2SqrTiled = permute(y2SqrTiled,[3 1 4 2 5]);
          y2SqrTiled = reshape(y2SqrTiled,[size(y2SqrTiled,1)*size(y2SqrTiled,2) size(y2SqrTiled,3)*size(y2SqrTiled,4) size(y2SqrTiled,5)]);
          
          v = imfilter(y2SqrTiled, params.muConvKernel );
          if params.muFeatureKernel
            vReduced = reshape(v,[size(v,1)*size(v,2) size(v,3)]);
            u = vReduced * params.muFeatureKernel;
          else
            u = v;
          end
          Psi_mu = (u + params.muEpsilon) .^ params.muPvalue;
          cost_mu = cost_mu + sum(Psi_mu(:));
        end
        
        if params.nu~=0
          y2block = reshape(y2, [size(y2,1)*size(y2,2) size(y2,3)]);
          d_EMA_y2_pairs = d_EMA_y2_pairs + (y2block' * y2block);
          d_EMA_y2 = d_EMA_y2 + sum(y2block,1);
          d_EMA_y2_square = d_EMA_y2_square + sum(y2block.^2,1);
        end
        
        if params.maskingNoiseFractionHidden
          tmp=randperm(numel(y2));
          numelToSetZero = round(params.maskingNoiseFractionHidden*numel(y2));
          y2Masked = y2;
          y2Masked(tmp(1:numelToSetZero)) = 0;
          if params.maskingNoiseFractionHiddenRescale
            realRatioSetZero = numelToSetZero/numel(y2Masked);
            y2Masked = y2Masked / (1-realRatioSetZero);
          end
          if params.maskingNoiseFractionHiddenEffectGradient
            y2 = y2Masked;
          end
        end
        
        if params.useTiledConv
          if params.maskingNoiseFractionHidden
            y2Tiled = reshape(y2Masked,[y2NumTilesX y2NumTilesY this.tiledConvF.tileSizeX this.tiledConvF.tileSizeY this.tiledConvF.fOut]);            
          else
            y2Tiled = reshape(y2,[y2NumTilesX y2NumTilesY this.tiledConvF.tileSizeX this.tiledConvF.tileSizeY this.tiledConvF.fOut]);
          end
          y2Tiled = permute(y2Tiled,[3 1 4 2 5]);
          y2Tiled = reshape(y2Tiled,[size(y2Tiled,1)*size(y2Tiled,2) size(y2Tiled,3)*size(y2Tiled,4) size(y2Tiled,5)]);
          
          z3 = this.tiledConvB.convW(reshape(W2,this.dimW2),y2Tiled,'valid');
          z3NumTilesX = size(z3,1);
          z3NumTilesY = size(z3,2);
          z3 = permute(z3,[3 1 4 2 5]);
          z3 = reshape(z3,[size(z3,1)*size(z3,2) size(z3,3)*size(z3,4) size(z3,5)]);
        else
          if params.maskingNoiseFractionHidden
            z3 = conv3d(W2,y2Masked);
          else
            z3 = conv3d(W2,y2);
          end
        end
        
        if ~params.noBiasBack
            z3 = bsxfun(@plus,z3,reshape(b2,[1 1 numel(b2)]));
        end
        if params.reconstrSigmoid
          a3 = Autoencoder.sigmoid(z3);
        elseif params.reconstrRectifiedLinear
          a3=max(z3,0);
        else
          a3 = z3;
        end
        
        if ~noOptout
          %calc stats per feature type:
          sumActivation_z2 = sumActivation_z2 + sum(sum(z2,1),2);
          sumActivationSquared_z2 = sumActivationSquared_z2 + sum(sum(z2.^2,1),2);
          sumActivationFourthMoment_z2 = sumActivationFourthMoment_z2 + sum(sum(z2.^4,1),2);
          numActive_z2 = numActive_z2 + sum(sum(z2>0,1),2);
          
          sumActivation_a2 = sumActivation_a2 + sum(sum(a2,1),2);
          sumActivationSquared_a2 = sumActivationSquared_a2 + sum(sum(a2.^2,1),2);
          sumActivationFourthMoment_a2 = sumActivationFourthMoment_a2 + sum(sum(a2.^4,1),2);
          numActive_a2 = numActive_a2 + sum(sum(a2>0,1),2);
          
          sumActivation_y2 = sumActivation_y2 + sum(sum(y2,1),2);
          sumActivationSquared_y2 = sumActivationSquared_y2 + sum(sum(y2.^2,1),2);
          sumActivationFourthMoment_y2 = sumActivationFourthMoment_y2 + sum(sum(y2.^4,1),2);
          numActive_y2 = numActive_y2 + sum(sum(y2>0,1),2);
          
          numPerFeatureType_layer2 = numPerFeatureType_layer2 + size(a2,1)*size(a2,2);
          
        end
        
        %now cut out the center of the input data to match layer a3
        cut1=(size(data,1)-size(a3,1))/2;
        cut2=(size(data,2)-size(a3,2))/2;
        
        if params.useCrossEntropyError
          costAutoenc = costAutoenc + sum(reshape( - data(cut1+1:end-cut1,cut2+1:end-cut2,:,k).*log(a3) - (1-data(cut1+1:end-cut1,cut2+1:end-cut2,:,k)).*log(1-a3) ,[1 numel(a3)]));
        elseif params.useLinearAutoencError
          costAutoenc = costAutoenc + sum(abs(reshape(a3 - data(cut1+1:end-cut1,cut2+1:end-cut2,:,k),[1 numel(a3)])));
        else %useSquaredError
          costAutoenc = costAutoenc + sum(reshape(a3 - data(cut1+1:end-cut1,cut2+1:end-cut2,:,k),[1 numel(a3)]).^2 / 2);
        end
        
        if nargout>1 && ~skipGrad
          
          if params.reconstrSigmoid && params.useCrossEntropyError
            %% special case, where it is easier to directly calculate the error of z3:
            delta3 = - (data(cut1+1:end-cut1,cut2+1:end-cut2,:,k) - a3) * params.alpha;
          else
            %% first calculate the reconstruction error:
            if params.useCrossEntropyError
              delta3 = ((data(cut1+1:end-cut1,cut2+1:end-cut2,:,k) - a3) ./ ( (a3-1).*a3 )) * params.alpha;
            elseif params.useLinearAutoencError
              delta3 = - sign(data(cut1+1:end-cut1,cut2+1:end-cut2,:,k) - a3) * params.alpha;
            else
              delta3 = - (data(cut1+1:end-cut1,cut2+1:end-cut2,:,k) - a3) * params.alpha;
            end

            %% now backpropagate the error through the activation function:
            if params.reconstrSigmoid
              delta3 = delta3 .* (a3.*(1-a3));
            elseif params.reconstrRectifiedLinear
              delta3 = delta3 .* sign(a3);
            else
              %do nothing
            end
          end
          
          delta3Full = zeros(size(delta3)+[(this.backwardRfSize(1)-1)*2 (this.backwardRfSize(2)-1)*2 0]);%zeros(size(data(:,:,:,k)));
          cut1Delta3Full = (this.backwardRfSize(1)-1);
          cut2Delta3Full = (this.backwardRfSize(2)-1);
          delta3Full(cut1Delta3Full+1:end-cut1Delta3Full,cut2Delta3Full+1:end-cut2Delta3Full,:) = delta3;
          
          delta2 = zeros(size(y2));
          
          if params.alpha~=0
            if params.useTiledConv
              
              W2inv = this.tiledConvB.invertW(W2);
              tmp = this.tiledConvBInv.convW(W2inv,delta3Full,'valid');
              tmp = reshape(tmp,[size(tmp,1) size(tmp,2) size(tmp,3)*size(tmp,4)*size(tmp,5)]);
              delta2 = bsxfun(@plus,delta2,tmp);
              
            else
              delta2 = bsxfun(@plus,delta2,conv3d(permute(flipdim(flipdim(W2,1),2),[1 2 4 3]),delta3Full));
            end
            if ~noOptout
              sumSquaredDelta2Autoenc = sumSquaredDelta2Autoenc + sum(delta2(:).^2);
              sumDelta2Autoenc = sumDelta2Autoenc + sum(delta2(:));

            end
          end
          
          if params.beta~=0
            delta2Sparse = params.beta*(- params.sparsityParam./avgAct2 + (1-params.sparsityParam)./(1-avgAct2) )  / numUnits2;
            delta2 = bsxfun(@plus,delta2,delta2Sparse);
            if ~noOptout
              sumSquaredDelta2Sparse = sumSquaredDelta2Sparse + sum(delta2Sparse(:).^2);
              sumDelta2Sparse = sumDelta2Sparse + sum(delta2Sparse(:));
            end
          end
          
          if params.nu~=0
            % dPsi_dy2 consists of three partial derivatives:
            
            % 1. dPsi/dEMA_y2_square * dEMA_y2_square/dy2
            delta2Decorr = bsxfun(@times, 2 * reshape(delta_EMA_y2_square,[1 1 size(y2,3)]), y2);
            
            % 2. dPsi/dEMA_y2 * dEMA_y2/dy2
            delta2Decorr = bsxfun(@plus, delta2Decorr, reshape(delta_EMA_y2,[1 1 size(y2,3)]));
            
            % 3. dPsi/delta_EMA_y2_pairs * delta_EMA_y2_pairs/dy2
            %delta2Decorr = delta2Decorr + 2 * reshape(sum( bsxfun(@times, reshape(delta_EMA_y2_pairs,[1 1 size(y2,3) size(y2,3)]), y2) ,3),[size(y2,1) size(y2,2) size(y2,3)]);
            % or the same more efficient:
            delta2Decorr = delta2Decorr + 2 * reshape( reshape(y2,[size(y2,1)*size(y2,2) size(y2,3)]) * delta_EMA_y2_pairs ,[size(y2,1) size(y2,2) size(y2,3)]);

            % all three terms above contain derivatives of the EMA. These
            % derivatives contain the EMAconst:
            delta2Decorr = delta2Decorr *  (1-params.nuEMAconst) / (size(y2,1)*size(y2,2));
            
            delta2  = bsxfun(@plus,delta2,delta2Decorr);
          end
          
          if params.mu~=0
            deltaU = params.muPvalue * ((u + params.muEpsilon).^(params.muPvalue-1));
            if params.muFeatureKernel
              deltaV = deltaU * params.muFeatureKernel;
              deltaV = reshape(deltaV, [size(v,1) size(v,2) size(v,3)]);
            else
              deltaV = deltaU;
            end
            deltaYSqr = imfilter(deltaV, params.muConvKernel, 'conv' );
            deltaYSqr = reshape(deltaYSqr,[this.tiledConvF.tileSizeX y2NumTilesX this.tiledConvF.tileSizeY y2NumTilesY this.tiledConvF.fOut]);
            deltaYSqr = permute(deltaYSqr,[2 4 1 3 5]);
            deltaYSqr = reshape(deltaYSqr,[y2NumTilesX y2NumTilesY this.tiledConvF.tileSizeX*this.tiledConvF.tileSizeY*this.tiledConvF.fOut]);
            delta2mu = 2 * deltaYSqr .* y2 / numUnits2;
            delta2 = delta2 + delta2mu;
            if ~noOptout
              sumSquaredDelta2mu = sumSquaredDelta2mu + sum(delta2mu(:).^2);
              sumDelta2mu = sumDelta2mu + sum(delta2mu(:));
            end
          end
          
          if params.gamma~=0
            deltaSparse = zeros(size(delta2,1)*size(delta2,2),size(delta2,3));
            if params.gammaPvalue==1
              for f2=1:size(z2,3)
                %TODO: change to im2col:
                dPsiToposparse_da2 = bsxfun(@times,y2reduced(:,f2),this.nghMat(f2,:)) ./ Psi_sparse;
                deltaSparse(:,f2) = deltaSparse(:,f2) + sum( dPsiToposparse_da2, 2 );
              end
            else
              exponent = (0.5 * params.gammaPvalue - 1)/(0.5 * params.gammaPvalue);
              tmpSparse = params.gammaPvalue * Psi_sparse.^exponent;
              for f2=1:size(z2,3)
                %TODO: change to im2col:
                dPsiToposparse_dy2 = bsxfun(@times,y2reduced(:,f2),this.nghMat(f2,:)) .* tmpSparse;
                deltaSparse(:,f2) = deltaSparse(:,f2) + sum( dPsiToposparse_dy2, 2 );
              end
            end
            delta2Toposparse = params.gamma*reshape(deltaSparse,size(delta2)) / numUnits2;
            delta2 = bsxfun(@plus,delta2, delta2Toposparse);
            if ~noOptout
              sumSquaredDelta2Toposparse = sumSquaredDelta2Toposparse + sum(delta2Toposparse(:).^2);
              sumDelta2Toposparse = sumDelta2Toposparse + sum(delta2Toposparse(:));
            end
          end
          
          % here is delta2 == dPsi_dy2 
          
          if params.useSoftmax
            % Here is delta2=dPsi_dy2 where y2=exp(a2)/sum(exp(a2),3)
            % where a2 is the activation after sigmoid or linear threshold
            tmpDelta2y2 = delta2 .* y2;
            
            % In several lines of code:
            %tmpDelta2y2SumWithoutOne = bsxfun(@minus, sum( tmpDelta2y2, 3), tmpDelta2a2);
            %tmpDelta2smax_notSameCoord = tmpDelta2y2SumWithoutOne .* y2;
            %tmpDelta2smax_sameCoord = delta2 .* (y2.*(1-y2));
            %delta2 = tmpDelta2smax_sameCoord - tmpDelta2smax_notSameCoord;
            
            % The same in one line of code:
            delta2 = delta2 .* (y2.*(1-y2)) - bsxfun(@minus, sum( tmpDelta2y2, 3), tmpDelta2y2) .* y2;
          elseif params.useNormMean
            
            delta2 = bsxfun(@rdivide, bsxfun(@minus, bsxfun(@times,delta2,localNormalizer), sum(delta2.*a2,3) ) , localNormalizer.^2);
            
          end
          
          % here is delta2 == dPsi_da2 
          
          if params.useRectifiedLinear
            delta2 = delta2 .* sign(a2);
          elseif params.useRectifiedLog
            delta2(a2~=0) = delta2(a2~=0) ./ (z2(a2~=0)+1);
            delta2(a2==0) = 0;
          elseif params.useNoActFcn
            %do nothing
          elseif params.useBinaryLinearInterp
             delta2(a2==0) = 0;
             delta2(a2==1) = 0;
             delta2 = bsxfun(@rdivide,delta2,reshape(EMA_z2pos,[1 1 size(delta2,3)])*params.binaryLinearInterpScaling);
          else
            delta2 = delta2 .* (a2.*(1-a2));
          end
          
          % here is delta2 == dPsi_dz2
          
          if params.useTiledConv
            % TODO: could be inlined to save memory:
            W1gradTmp = this.tiledConvF.backpropErrorW(curData,reshape(delta2,[y2NumTilesX y2NumTilesY this.tiledConvF.tileSizeX this.tiledConvF.tileSizeY this.tiledConvF.fOut]));
            W1grad = W1grad + reshape(W1gradTmp,[size(W1gradTmp,1) size(W1gradTmp,2) size(W1gradTmp,3) size(W1gradTmp,4)*size(W1gradTmp,5)*size(W1gradTmp,6)]);
          else
            W1grad = W1grad + Autoencoder.conv3d3(delta2,curData);
          end
          
          if params.useTiledConv
            delta3reshaped = permute(reshape(delta3,[size(delta3,1)/z3NumTilesX z3NumTilesX size(delta3,2)/z3NumTilesY z3NumTilesY size(delta3,3)]),[2 4 1 3 5]);
%             reshape(y2,[size(y2,1) size(y2,2) y2NumTilesX y2NumTilesY this.tiledConvF.fOut])
            deltaW2 = this.tiledConvB.backpropErrorW(y2Tiled,delta3reshaped);
            deltaW2 = reshape(deltaW2,[size(deltaW2,1) size(deltaW2,2) size(deltaW2,3) size(deltaW2,4)*size(deltaW2,5)*size(deltaW2,6)]);
          else
            deltaW2 = Autoencoder.conv3d3(delta3,y2);
          end
          
          if params.fixL2WeightBack
            if params.fixL2WeightBackPoolPerHiddenNeuron
              deltaW2tilde = bsxfun(@rdivide, deltaW2, W2tildeNormalizer) - bsxfun(@times, W2tilde, sum(sum(sum( deltaW2 .* W2tilde, 1),2),4 ) ./ W2tildeNormalizer.^3 );
            else
              deltaW2tilde = bsxfun(@rdivide, deltaW2, W2tildeNormalizer) - bsxfun(@times, W2tilde, sum(sum(sum( deltaW2 .* W2tilde, 1),2),3 ) ./ W2tildeNormalizer.^3 );
            end
            W2grad = W2grad + deltaW2tilde;
          else
            W2grad = W2grad + deltaW2;
          end
          
          b1grad = b1grad + squeeze(sum(sum(delta2,2),1));
          if ~params.noBiasBack
              b2grad = b2grad + squeeze(sum(sum(delta3,2),1));
          end
          
        end
      end
      
      if params.nu~=0
        d_EMA_y2_pairs = d_EMA_y2_pairs / (size(y2,1)*size(y2,2)*m);
        d_EMA_y2 = d_EMA_y2 / (size(y2,1)*size(y2,2)*m);
        d_EMA_y2_square = d_EMA_y2_square / (size(y2,1)*size(y2,2)*m);
          
        EMA_y2_pairs = EMA_y2_pairs*params.nuEMAconst + d_EMA_y2_pairs *(1-params.nuEMAconst);
        EMA_y2 = EMA_y2*params.nuEMAconst + d_EMA_y2 * (1-params.nuEMAconst);
        EMA_y2_square = EMA_y2_square*params.nuEMAconst + d_EMA_y2_square * (1-params.nuEMAconst);
        
        %% again as in the beginning but now just for calculating the cost with updated values:
        % TODO: this could be optimized because the function getVariances
        % is called here and in the beginning of the next call to
        % autoencoderCost()
        [cov_y2, var_y2, var_pair] = this.getVariances();
        if ~isempty(params.nuCovWeightMat)
          cost_nu = params.nu * sum( params.nuCovWeightMat(:) .* cov_y2(:).^2 ./ var_pair(:) );
        else
          cost_nu = params.nu * sum( cov_y2(:).^2 ./ var_pair(:) );
        end
      end
      
      
      
      if params.gamma~=0
        avgPsi_sparse = avgPsi_sparse/(m*numUnits2);
        costTopoSparse = params.gamma*sum(avgPsi_sparse);
      else
        costTopoSparse = 0;
      end
      if params.mu~=0
        cost_mu = params.mu*cost_mu/(m*numUnits2);
      end
      costAutoenc = params.alpha*costAutoenc / m;
      b1grad = b1grad / m;
      b2grad = b2grad / m;
      W1grad = W1grad / m;
      W2grad = W2grad / m;
      if ~noOptout
        optout.sumSquaredW1gradAutoencSparse = sum(W1grad(:).^2);
        optout.sumSquaredW2gradAutoencSparse = sum(W2grad(:).^2);
      end
      
      if params.lambda~=0
        if ~noOptout
          optout.sumSquaredW1gradRegL2 = params.lambda * sum(W1(:).^2);
          optout.sumSquaredW2gradRegL2 = (params.lambda * params.lambdaBackScale) * sum(W2(:).^2);
        end
        W1grad = W1grad + params.lambda * W1;
        W2grad = W2grad + (params.lambda * params.lambdaBackScale) * W2;
      end
      
      if params.lambdaL1~=0
        if params.smoothLambdaL1WithEpsilon>0
          W1gradRegL1 = params.lambdaL1 * W1./sqrt(params.smoothLambdaL1WithEpsilon+W1.^2);
          W2gradRegL1 = (params.lambdaL1*params.lambdaL1BackScale) * W2./sqrt(params.smoothLambdaL1WithEpsilon+W2.^2);
        else
          W1gradRegL1 = params.lambdaL1 * sign(W1);
          W2gradRegL1 = (params.lambdaL1*params.lambdaL1BackScale) * sign(W2);
        end
        if params.scaleLambdaL1WithDistExponent
          W1gradRegL1 = bsxfun(@times,distMatW1,W1gradRegL1);
          W2gradRegL1 = bsxfun(@times,distMatW2,W2gradRegL1);
        end
        W1grad = W1grad + W1gradRegL1;
        W2grad = W2grad + W2gradRegL1;
        if ~noOptout
          optout.sumSquaredW1gradRegL1 = sum(W1gradRegL1(:).^2);
          optout.sumSquaredW2gradRegL1 = sum(W2gradRegL1(:).^2);
        end
      end
      
      if nargout > 1
        grad = [W1grad(:) ; W2grad(:) ; b1grad(:) ; b2grad(:)];
      end
      cost = costAutoenc + costRegularize + costRegularizeL1 + costSparseness + costTopoSparse + cost_mu + cost_nu;
      
      if ~noOptout
        optout.costAutoenc = costAutoenc;
        optout.costRegularize = costRegularize;
        optout.costRegularizeL1 = costRegularizeL1;
        optout.costSparseness = costSparseness;
        optout.cost_mu = cost_mu;
        optout.costTopoSparse = costTopoSparse;
        %   optout.grad.b1grad = b1grad;
        %   optout.grad.b2grad = b2grad;
        %   optout.grad.W1grad = W1grad;
        %   optout.grad.W2grad = W2grad;
        
        optout.hardSparseness_z2 = numActive_z2(:)./numPerFeatureType_layer2(:);
        optout.meanActivation_z2 = sumActivation_z2(:)./numPerFeatureType_layer2(:);
        optout.stdActivation_z2 = sqrt( sumActivationSquared_z2(:)./numPerFeatureType_layer2(:) - (sumActivation_z2(:)./numPerFeatureType_layer2(:)).^2 );
        optout.hardSparseness_a2 = numActive_a2(:)./numPerFeatureType_layer2(:);
        optout.meanActivation_a2 = sumActivation_a2(:)./numPerFeatureType_layer2(:);
        optout.stdActivation_a2 = sqrt( sumActivationSquared_a2(:)./numPerFeatureType_layer2(:) - (sumActivation_a2(:)./numPerFeatureType_layer2(:)).^2 );
        optout.hardSparseness_y2 = numActive_y2(:)./numPerFeatureType_layer2(:);
        optout.meanActivation_y2 = sumActivation_y2(:)./numPerFeatureType_layer2(:);
        optout.stdActivation_y2 = sqrt( sumActivationSquared_y2(:)./numPerFeatureType_layer2(:) - (sumActivation_y2(:)./numPerFeatureType_layer2(:)).^2 );
        optout.kurtosisActivation_y2 = (sumActivationFourthMoment_y2(:)/numPerFeatureType_layer2(:)) ./ (sumActivationSquared_y2(:)/numPerFeatureType_layer2(:)).^2 - 3;

        optout.rfBrightness = mean(reshape(W1,[size(W1,1)*size(W1,2)*size(W1,3) size(W1,4)]),1);
        
        optout.meanW1 = mean(W1(:));
        optout.stdW1 = std(W1(:));
        optout.meanW2 = mean(W2(:));
        optout.stdW2 = std(W2(:));
        
        optout.meanb1 = mean(b1(:));
        optout.stdb1 = std(b1(:));
        optout.meanb2 = mean(b2(:));
        optout.stdb2 = std(b2(:));
        
        if ~skipGrad
          if params.alpha~=0
            optout.meanSquaredDelta2Autoenc = sumSquaredDelta2Autoenc / (m*numel(delta2));
            optout.meanDelta2Autoenc = sumDelta2Autoenc / (m*numel(delta2));
            optout.stdDelta2Autoenc = optout.meanSquaredDelta2Autoenc - optout.meanDelta2Autoenc^2;
          end
          if params.beta~=0
            optout.meanSquaredDelta2Sparse = sumSquaredDelta2Sparse / (m*numel(delta2));
            optout.meanDelta2Sparse = sumDelta2Sparse / (m*numel(delta2));
            optout.stdDelta2Sparse = optout.meanSquaredDelta2Sparse - optout.meanDelta2Sparse^2;
          end
          if params.gamma~=0
            optout.meanSquaredDelta2Toposparse = sumSquaredDelta2Toposparse / (m*numel(delta2));
            optout.meanDelta2Toposparse = sumDelta2Toposparse / (m*numel(delta2));
            optout.stdDelta2Toposparse = optout.meanSquaredDelta2Toposparse - optout.meanDelta2Toposparse^2;
          end
          if params.mu~=0
            optout.meanSquaredDelta2mu = sumSquaredDelta2mu / (m*numel(delta2));
            optout.meanDelta2mu = sumDelta2mu / (m*numel(delta2));
            optout.stdDelta2mu = optout.meanSquaredDelta2mu - optout.meanDelta2mu^2;
          end
        end
      end
      
      if noStateChange
        minibatchIds = temp_minibatchIds;
        EMA_z2pos =  temp_EMA_z2pos;
        EMA_y2_pairs =  temp_EMA_y2_pairs;
        EMA_y2_square =  temp_EMA_y2_square;
        EMA_y2 =  temp_EMA_y2;
      end

    end
    
    function [curData,z2,a2,y2,localNormalizer] = forwardPass(this,curData,W1,b1)
      global EMA_z2pos;
      
      params = this.params.Autoencoder;
      if params.maskingNoiseFraction
        tmp=randperm(numel(curData));
        numelToSetZero = round(params.maskingNoiseFraction*numel(curData));
        curData(tmp(1:numelToSetZero)) = 0;
        if params.maskingNoiseFractionRescale
          realRatioSetZero = numelToSetZero/numel(curData);
          curData = curData / (1-realRatioSetZero);
        end
      end
      if params.gaussianNoiseSigmaInput
        if params.gaussianNoiseCorrelated
          if params.gaussianNoiseCorrelatedUniform
            curData = curData + params.gaussianNoiseSigmaInput * randn(size(curData)) * rand(1);
          else
            curData = curData + params.gaussianNoiseSigmaInput * randn(size(curData)) * randn(1);
          end
        else
          curData = curData + params.gaussianNoiseSigmaInput * randn(size(curData));
        end
      end
      
      if params.useTiledConv
        z2 = this.tiledConvF.convW(reshape(W1,this.dimW1),curData,'valid');
        z2 = reshape(z2,[size(z2,1) size(z2,2) size(z2,3)*size(z2,4)*size(z2,5)]);
      else
        z2 = conv3d(W1,curData);
      end
      z2=bsxfun(@plus,z2,reshape(b1,[1 1 numel(b1)]));
      if params.useRectifiedLinear
        a2=max(z2,0);
      elseif params.useRectifiedLog
        a2 = log(max(z2,0)+1);
      elseif params.useNoActFcn
        a2 = z2;
      elseif params.useBinaryLinearInterp
        a2 = min(max(bsxfun(@rdivide,z2,reshape(EMA_z2pos,[1 1 numel(EMA_z2pos)])*params.binaryLinearInterpScaling),0),1);
      else
        a2 = Autoencoder.sigmoid(z2);
      end
      
      if params.hiddenLinearBinaryThresholdUnitsFirst
        a2 = (a2>0);
      elseif params.hiddenSigmoidBinaryThresholdUnitsFirst
        a2 = (a2>0.5);
      end
      
      if params.useSoftmax
        localNormalizer = [];
        y2 = Autoencoder.softmax(a2);
      elseif params.useNormMean
        localNormalizer = sum(a2,3)+params.normMeanEpsilon;
        y2 = bsxfun(@rdivide,a2,localNormalizer);
      else
        localNormalizer = [];
        y2 = a2;
      end
      
      if params.hiddenLinearBinaryThresholdUnits
        y2 = (y2>0);
      elseif params.hiddenSigmoidBinaryThresholdUnits
        y2 = (y2>0.5);
      end
    end
    
    function [cov_y2, var_y2, var_pair] = getVariances(this)
      global EMA_y2_pairs;
      global EMA_y2_square;
      global EMA_y2;
      
      cov_y2 = (EMA_y2_pairs - EMA_y2' * EMA_y2);
      var_y2 = (EMA_y2_square - EMA_y2.^2);
      var_pair = (var_y2' * var_y2);
      
      %set diagonal to zero because we don't want to use them later when we sum over rows or columns
      cov_y2(logical(eye(size(cov_y2)))) = 0; 
      
      %exact value in the diagonal doesn't matter. it should also work with any postive value.. 
      %it is just used such that a division in the diagonal elements will result in 0:
      var_pair(logical(eye(size(var_pair)))) = 1; 
      
    end
  end
  
  methods (Static)
    function z2 = conv3d3(W,data)
      z2 = zeros([size(data,1)-size(W,1)+1 size(data,2)-size(W,2)+1 size(data,3) size(W,3)]);
      for f1=1:size(data,3)
        blockImg = im2col(data(:,:,f1),[size(W,1) size(W,2)]);
        blockFilters = reshape(W(:,:,:),[size(W,1)*size(W,2) size(W,3)]);
        z2(:,:,f1,:) = reshape(blockImg'*blockFilters,[size(z2,1) size(z2,2) 1 size(W,3)]);
      end
    end
    
    function sigm = sigmoid(x)
      sigm = 1 ./ (1 + exp(-x));
    end
    
    function smax = softmax(x)
      smax = bsxfun(@rdivide,exp(x),sum(exp(x),3));
    end
   
  end
  
  methods
    function plotLastIters( this )
      
      %       if this.numJobs > 1
      for i=1:this.numJobs
        tmp = dir(fullfile(num2str(i),'forwConnIter*.mat'));
        tmp = sort(cellfun(@(x) str2double(x(13:end-4)),{tmp.name}));
        lastIter=tmp(end);
        forwConn = load(fullfile(this.temppath,num2str(i),['forwConnIter' num2str(lastIter) '.mat']));
        sizeW = size(forwConn.W);
        forwConn.W = reshape(forwConn.W,[sizeW(1:3) prod(sizeW(4:6))]);
        sizeW = size(forwConn.W);
        plotdim = ceil(sqrt(sizeW(4)));
        forwConn.W = cat(4,forwConn.W,zeros([sizeW(1:3) plotdim.^2-sizeW(4)]));
        plotColorFeatures(  reshape(forwConn.W,[sizeW(1:3) plotdim plotdim]), true, fullfile(this.temppath,['jobid' num2str(i) 'forwConnIter' num2str(lastIter) '.png']), true );
      end
      %       else
      %         weights=load(fullfile(this.temppath,'currIter.mat'),'W1');
      %         if exist('whitening','var')
      %           weights.W1 = addWhiteningWeights(weights.W1, whitening.W);
      %         end
      %         savepath = fullfile(this.resultpath, 'weights.png');
      %         plotColorFeatures(  reshape(weights.W1,[12 12 3 10 size(weights.W1,4)/10]), true, savepath, true );
      %       end
      
    end
    
    function plotLastValidationStatistics( this, varname )
      if nargin<2
        disp(fieldnames(this.log.validation(end).iterationType.optout))
        varname = 'y2';
      end
      [B,IX] = sort(this.log.validation(end).iterationType.optout.(['meanActivation_' varname]));
     
      figure;
      subplot(2,1,1)
      errorbar(this.log.validation(end).iterationType.optout.(['meanActivation_' varname])(IX),this.log.validation(end).iterationType.optout.(['stdActivation_' varname])(IX))
      xlim([0 length(IX)+1])
      ylabel('mean activation')
      
      subplot(2,1,2)
      plot(this.log.validation(end).iterationType.optout.(['hardSparseness_' varname])(IX))
      xlim([0 length(IX)+1])
      xlabel('neuron id')
      ylabel('fraction of time active')
     
      figure;
      subplot(2,1,1)
      hist(this.log.validation(end).iterationType.optout.(['meanActivation_' varname]));
      xlabel('mean activation')
      ylabel('number of neurons')
      
      subplot(2,1,2)
      hist(this.log.validation(end).iterationType.optout.(['hardSparseness_' varname]));
      xlabel('fraction of time active')
      ylabel('number of neurons')
      
      
    end
    
    function plotValidationStatistics( this )
      fnames = fieldnames(this.log.validation(end).iterationType.optout);
      [Selection,ok] = listdlg('ListString',fnames);
      
      for t=1:length(this.log.validation)
        tmp = this.log.validation(t).iterationType.optout.(fnames{Selection});
        cumStat{t} = reshape(tmp,[1 numel(tmp)]);
        i(t) = this.log.validation(t).i;
      end
      cumStat = cell2mat(cumStat');
      
      figure;
      if size(cumStat,2)>1
        imagesc(i,1:size(cumStat,2),cumStat')
        title(fnames{Selection})
      else
        plot(i,cumStat)
        ylabel(fnames{Selection})
      end
      xlabel('iteration')
      
    end
    
    function plotIterationInfos( this, doValidationSet, fnames, jobids, plotPerJob )
      if ~isdir(this.resultpath)
        mkdir(this.resultpath);
      end
      
      if nargin<3
        fnames = [];
      end
      if nargin<4
        jobids = [];
      end
      if nargin<5
        plotPerJob = true;
      end
      
      if doValidationSet
        matchfilenames = 'validInfoIter.*.mat';
      else
        matchfilenames = 'trainInfoIter.*.mat';
      end
      
      if isempty(jobids)
        jobids = 1:this.numJobs;
      end
      
      
      optout = cell(1,this.numJobs);
      for itmp=1:length(jobids)
        i=jobids(itmp);
        if this.numJobs>1
          files = dir(fullfile(this.temppath,num2str(i)));
        else
          files = dir(fullfile(this.temppath));
        end
        counter = 1;
        iterationTypeAll = cell(size(files));
        for j=1:length(files)
          if ~isempty(regexp(files(j).name, matchfilenames, 'once'))
            if this.numJobs > 1
              iterationType = load(fullfile(this.temppath,num2str(i),files(j).name),'iterationType');
            else
              iterationType = load(fullfile(this.temppath,files(j).name),'iterationType');
            end
            if strcmp(iterationType.iterationType,'init')% || isempty(iterationType.iterationType.cost)
              continue;
            end
            iterationTypeAll{counter} = iterationType.iterationType;
            counter = counter + 1;
          end
        end
        iterationTypeAll(counter:end) = [];
        
        optout{i}.iteration = cellfun(@(x) x.iteration, iterationTypeAll);
        [optout{i}.iteration,sortIds] = sort(optout{i}.iteration);
        
        optout{i}.cost = cellfun(@(x) x.cost, iterationTypeAll);
        optout{i}.cost = optout{i}.cost(sortIds);
        
        optoutFields = fieldnames(iterationTypeAll{1}.optout);
        for k=1:length(optoutFields)
          optout{i}.(optoutFields{k}) = cellfun(@(x) x.optout.(optoutFields{k}), iterationTypeAll);
          optout{i}.(optoutFields{k}) = optout{i}.(optoutFields{k})(sortIds);
        end
        
      end
      
      optoutFields = fieldnames(optout{i});
      
      disp(optoutFields);
      if isempty(fnames)
        fnames = optoutFields;
      end
      
      
      if plotPerJob
        cmap = colormap('lines');
        for itmp=1:length(jobids)
          i=jobids(itmp);
          figure; clf;
          legendText = cell(0);
          for j=1:length(fnames)
            if isfield(optout{i},fnames{j})
              plot(optout{i}.iteration, optout{i}.(fnames{j}),'Color',cmap(j,:));
              hold on;
              legendText{end+1} = fnames{j}; %#ok<AGROW>
            end
          end
          hold off;
          xlabel('iteration')
          jobName = '';
          for k=1:size(this.paramComb,1)
            jobName = [jobName ' ' this.variableParams{k}{2} '=' num2str(this.paramComb{k,i})];
          end
          ylabel(jobName)
          legend(legendText);
        end
      else
        cmap = colormap('lines');
        for j=1:length(fnames)
          figure; clf;
          legendText = cell(0);
          for itmp=1:length(jobids)
            i=jobids(itmp);
            if isfield(optout{i},fnames{j})
              plot(optout{i}.iteration, optout{i}.(fnames{j}),'Color',cmap(i,:));
              hold on;
              legendText{end+1} = ''; %#ok<AGROW>
              for k=1:size(this.paramComb,1)
                legendText{end} = [legendText{end} ' ' this.variableParams{k}{2} '=' num2str(this.paramComb{k,i})];
              end
            end
          end
          hold off;
          xlabel('iteration')
          ylabel(fnames{j})
          legend(legendText);
        end
        
      end
      
    end
    
    function weights = addWhiteningWeights(weights, whiteningWeights)
      %TODO
      %origSize = size(weights);
      weights = reshape(weights,[]);
      weights = whiteningWeights * weights;
    end
    
    
    function [diff,numGrad,grad,diffs] = checkGradient(this, patches)
      
      %this.forwardRfSize = [this.params.Autoencoder.patchDimForward this.params.Autoencoder.patchDimForward size(patches,3)];
      %this.backwardRfSize = [this.params.Autoencoder.patchDimBackward this.params.Autoencoder.patchDimBackward this.params.Autoencoder.hiddenSize];
      
      this = this.setPatchsizes();
      
%       W1 = randn([this.forwardRfSize this.backwardRfSize(3)]);
%       W2 = randn([this.backwardRfSize this.forwardRfSize(3)]);
%       b1 = randn(this.backwardRfSize(3),1);
%       b2 = randn(this.forwardRfSize(3),1);
      W1 = randn(this.dimW1);
      W2 = randn(this.dimW2);
      b1 = randn(this.dimW1(4:6));
      b2 = randn([ 1 1 this.dimW2(6)]);
      theta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
      
      if this.params.Autoencoder.gamma~=0
        if isempty(this.params.Autoencoder.topoNghFcn)
          this.nghMat = linearDecoderNghMat( this.params.Autoencoder.topoNgh,this.params.Autoencoder.topoGridDims,this.params.Autoencoder.topoPeriodicBoundary);
        else
          this.nghMat = linearDecoderNghWinMat( this.params.Autoencoder.topoNgh,this.params.Autoencoder.topoGridDims,this.params.Autoencoder.topoPeriodicBoundary, this.params.Autoencoder.topoNghFcn);
        end
      end
      
      if this.params.Autoencoder.nu
        % initialize EMA of activitives:
        for k=1:100
          cost = this.autoencoderCost(theta, patches);
        end
      end
      
      globStream = Autoencoder.getGlobRng();
      globStreamInitState = globStream.State;
      [~, grad] = this.autoencoderCost(theta, patches, false, true, true);
      
      globStream.State = globStreamInitState;
      numGrad = Autoencoder.computeNumericalGradient(@(p) autoencoderCost(this, p, patches, true, true, true), theta);
      
      diff = norm(numGrad-grad)/norm(numGrad+grad);
      disp(['normalized diff grad: ' num2str(diff)]); 
      
      %% calc diff per var
      [grads.W1, grads.W2, grads.b1, grads.b2] = this.decodeTheta(grad);
      [numGrads.W1, numGrads.W2, numGrads.b1, numGrads.b2] = this.decodeTheta(numGrad);
      
      diffs.b1 = norm(numGrads.b1(:)-grads.b1(:))/norm(numGrads.b1(:)+grads.b1(:));
      diffs.b2 = norm(numGrads.b2(:)-grads.b2(:))/norm(numGrads.b2(:)+grads.b2(:));
      diffs.W1 = norm(numGrads.W1(:)-grads.W1(:))/norm(numGrads.W1(:)+grads.W1(:));
      diffs.W2 = norm(numGrads.W2(:)-grads.W2(:))/norm(numGrads.W2(:)+grads.W2(:));
      
      disp(['normalized diffs.b1=' num2str(diffs.b1)]); 
      disp(['normalized diffs.b2=' num2str(diffs.b2)]); 
      disp(['normalized diffs.W1=' num2str(diffs.W1)]); 
      disp(['normalized diffs.W2=' num2str(diffs.W2)]); 
    end
  end
  
  methods (Static)
    
    
    function unittest()
      
      this = Autoencoder();
      this.params.Autoencoder.patchDimForward = 5;
      this.params.Autoencoder.patchDimBackward = 3;
      this.params.Autoencoder.hiddenSize = 16;
      this.params.Autoencoder.useSoftmax = true;
      this.params.Autoencoder.useRectifiedLinear = false;
      this.params.Autoencoder.useNoActFcn = true;
      this.params.Autoencoder.sparsityParam = 0.035;
      this.params.Autoencoder.alpha = 1;
      this.params.Autoencoder.lambda = 0;
      this.params.Autoencoder.lambdaL1 = 0;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0.1;
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 2;
      this.params.Autoencoder.beta = 0;
      this.params.Autoencoder.gamma = 1;
      this.params.Autoencoder.gammaPvalue = 0.1;
      this.params.Autoencoder.maskingNoiseFraction = 0.25;
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0.5;
      this.params.Autoencoder.topoNghFcn = [];%@(x) gausswin(x,3);
      this.params.Autoencoder.topoNgh = [3 3];
      this.params.Autoencoder.topoGridDims = [4 4];
      this.params.Autoencoder.topoPeriodicBoundary = [true true];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      patches = randn([20 20 3 3]);
      this.checkGradient(patches);
      
    end
    
    
    
    function [diff,numGrad,grad,diffs] = unittestDecorr()
      
      this = Autoencoder();
      
      this.params.Autoencoder.patchDimForward = 3;
      this.params.Autoencoder.patchDimBackward = 3;
      this.params.Autoencoder.hiddenSize = 10;
      this.params.Autoencoder.useTiledConv = false; %if true: tileSize=patchDimForward-1, fOut=hiddenSize      
      this.params.Autoencoder.inputSubsampling = 1;
      
      % PROCESSING STAGES:
      this.params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
      this.params.Autoencoder.maskingNoiseFractionRescale = true;        % if inputs are rescaled to the variance before masking
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      this.params.Autoencoder.gaussianNoiseCorrelated = false;
      this.params.Autoencoder.gaussianNoiseCorrelatedUniform = false;

      this.params.Autoencoder.useRectifiedLinear = false; %instead of sigmoid
      this.params.Autoencoder.useRectifiedLog = false; %instead of sigmoid
      this.params.Autoencoder.useNoActFcn = false; %instead of sigmoid
      this.params.Autoencoder.useBinaryLinearInterp = false;
      this.params.Autoencoder.binaryLinearInterpScaling = 1;
      this.params.Autoencoder.binaryLinearInterpEMAconst = 0.99;
      
      this.params.Autoencoder.hiddenLinearBinaryThresholdUnitsFirst = false; %these are not considered when calculating gradients
      this.params.Autoencoder.hiddenSigmoidBinaryThresholdUnitsFirst = false; %these are not considered when calculating gradients

      this.params.Autoencoder.useSoftmax = false;
      this.params.Autoencoder.useNormMean = false;
      this.params.Autoencoder.normMeanEpsilon = 1e-6;
      
      this.params.Autoencoder.hiddenLinearBinaryThresholdUnits = false; %these are not considered when calculating gradients
      this.params.Autoencoder.hiddenSigmoidBinaryThresholdUnits = false; %these are not considered when calculating gradients
     
      this.params.Autoencoder.maskingNoiseFractionHidden = 0;        % fraction of inputs to set to 0
      this.params.Autoencoder.maskingNoiseFractionHiddenRescale = false;        % if inputs are rescaled to the variance before masking
      this.params.Autoencoder.maskingNoiseFractionHiddenEffectGradient = true; % if the masked y2 should be used when backpropagating errors
      
      this.params.Autoencoder.noBiasBack = false;

      this.params.Autoencoder.reconstrSigmoid = false;
      this.params.Autoencoder.reconstrRectifiedLinear = false;

      % OPTIMIZATION PARAMS:
      this.params.Autoencoder.alpha = 0.1;        % weight of Autoencoder 
      this.params.Autoencoder.useLinearAutoencError = false;
      this.params.Autoencoder.useCrossEntropyError = false;
      this.params.Autoencoder.beta = 0;         % weight of sparsity penalty term
      this.params.Autoencoder.sparsityParam = 0.035; % desired average activation of the hidden units.
      this.params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
      this.params.Autoencoder.gammaPvalue = 1;        % Psi_gamma = gamma * sqrt( NghMat * y2.^2 + epsilon ) ^ gammaPvalue
      this.params.Autoencoder.topoNghFcn = [];
      this.params.Autoencoder.topoNgh = [3 3];
      this.params.Autoencoder.topoGridDims = [20 20];
      this.params.Autoencoder.topoPeriodicBoundary = [true true];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      this.params.Autoencoder.mu = 0;        % weight of modular topographic sparsity penalty term
      this.params.Autoencoder.muPvalue = 1;        % Psi_mu = mu * ( muFeatureKernel * (muConvKernel * y2.^2) + muEpsilon ) ^ muPvalue
      this.params.Autoencoder.muConvKernel = ones(3,3);
      this.params.Autoencoder.muFeatureKernel = []; % if empty then no combination in feature dimension, otherwise it should be a matrix of size (hiddenSize,hiddenSize)
      this.params.Autoencoder.muEpsilon = 1e-2;
      this.params.Autoencoder.nu = 1;
      this.params.Autoencoder.nuEMAconst = 0.9;
      this.params.Autoencoder.nuCovWeightMat = diag(ones(1,9),1)+diag(ones(1,9),-1)+diag(ones(1,10),0);
      this.params.Autoencoder.lambda = 0;    % weight of L2-decay parameter
      this.params.Autoencoder.lambdaBackScale = 1;    % weight of L2-decay parameter scaling for backward conenctions
      this.params.Autoencoder.lambdaL1 = 0 ;    % weight of L1-decay parameter
      this.params.Autoencoder.lambdaL1BackScale = 1;    % weight of L1-decay parameter
      this.params.Autoencoder.fixL2WeightBack = false; %if true, remember to also set both backScale to 0
      this.params.Autoencoder.fixL2WeightBackPoolPerHiddenNeuron = false;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0 ; %if not zero then smooth it
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 0 ;
    
      patches = randn([20 20 3 2]);
      
      [diff,numGrad,grad,diffs] = this.checkGradient(patches);
      
    end
    
    function [diff,numGrad,grad] = unittestReLU()
      
      this = Autoencoder();
      this.params.Autoencoder.patchDimForward = 5;
      this.params.Autoencoder.patchDimBackward = 3;
      this.params.Autoencoder.hiddenSize = 16;
      this.params.Autoencoder.useSoftmax = false;
      this.params.Autoencoder.useNormMean = true;
      this.params.Autoencoder.useRectifiedLinear = true;
      this.params.Autoencoder.useNoActFcn = false;
      this.params.Autoencoder.sparsityParam = 0.035;
      this.params.Autoencoder.alpha = 1;
      this.params.Autoencoder.lambda = 0.1;
      this.params.Autoencoder.lambdaL1 = 0;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0.1;
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 2;
      this.params.Autoencoder.beta = 0;
      this.params.Autoencoder.maskingNoiseFraction = 0.25;
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      
      this.params.Autoencoder.topoNghFcn = [];
      this.params.Autoencoder.topoNgh = [1];
      this.params.Autoencoder.topoGridDims = [16];
      this.params.Autoencoder.topoPeriodicBoundary = [false];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      this.params.Autoencoder.gamma = 1;
      this.params.Autoencoder.gammaPvalue = 0.1;
      
      this.params.Autoencoder.useLinearAutoencError = false;
      patches = randn([20 20 3 3]);
      [diff,numGrad,grad,diffs] = this.checkGradient(patches);
      
    end
    
    
    function unittestSigmoidReconstr()
      
      this = Autoencoder();
      
      this.params.Autoencoder.patchDimForward = 5;
      this.params.Autoencoder.patchDimBackward = 5;
      this.params.Autoencoder.hiddenSize = 16;
      
      this.params.Autoencoder.useRectifiedLinear = false;
      this.params.Autoencoder.useNoActFcn = false;
      
      this.params.Autoencoder.useSoftmax = false;
      this.params.Autoencoder.useNormMean = false;
      
      this.params.Autoencoder.reconstrSigmoid = true;
      this.params.Autoencoder.reconstrRectifiedLinear = false;

      this.params.Autoencoder.useLinearAutoencError = false;
      this.params.Autoencoder.useCrossEntropyError = true;
      
      this.params.Autoencoder.sparsityParam = 0.035;
      this.params.Autoencoder.alpha = 1;
      this.params.Autoencoder.lambda = 1;
      this.params.Autoencoder.lambdaL1 = 1;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0.1;
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 2;
      this.params.Autoencoder.beta = 1;
      this.params.Autoencoder.gamma = 1;
      this.params.Autoencoder.gammaPvalue = 0.1;
      this.params.Autoencoder.maskingNoiseFraction = 0;
      this.params.Autoencoder.maskingNoiseFractionHidden = 0.1;
      this.params.Autoencoder.maskingNoiseFractionHiddenRescale = false;

      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      this.params.Autoencoder.topoNghFcn = @(x) gausswin(x,3);
      this.params.Autoencoder.topoNgh = [3 3];
      this.params.Autoencoder.topoGridDims = [4 4];
      this.params.Autoencoder.topoPeriodicBoundary = [true true];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      
      
      patches = rand([20 20 3 3]);
      this.checkGradient(patches);
      
    end
    
    function [diff,numGrad,grad,diffs] = unittestTiledConv()
      
      this = Autoencoder();
      
      this.params.Autoencoder.inSamplesDims = [18 18 3 3];
      patches = rand(this.params.Autoencoder.inSamplesDims);
      
      this.params.Autoencoder.patchDimForward = 4;
      this.params.Autoencoder.patchDimBackward = 4; %tileSize=patchDim-1
      this.params.Autoencoder.hiddenSize = 2;
      this.params.Autoencoder.useTiledConv = true;
      
      this.params.Autoencoder.useRectifiedLinear = true;
      this.params.Autoencoder.useNoActFcn = false;
      
      this.params.Autoencoder.useSoftmax = false;
      this.params.Autoencoder.useNormMean = false;
      
      this.params.Autoencoder.reconstrSigmoid = false;
      this.params.Autoencoder.reconstrRectifiedLinear = false;

      this.params.Autoencoder.useLinearAutoencError = false;
      this.params.Autoencoder.useCrossEntropyError = false;
      
      this.params.Autoencoder.sparsityParam = 0.035;
      this.params.Autoencoder.alpha = 0.01;
      this.params.Autoencoder.lambda = 0.01;
      this.params.Autoencoder.lambdaL1 = 0;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0.1;
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 2;
      this.params.Autoencoder.beta = 0;
      this.params.Autoencoder.gamma = 0;
      this.params.Autoencoder.gammaPvalue = 0.1;
      this.params.Autoencoder.mu = 1;        % weight of modular topographic sparsity penalty term
      this.params.Autoencoder.muPvalue = 0.2;        % Psi_mu = mu * ( muFeatureKernel * (muConvKernel * y2.^2) + muEpsilon ) ^ muPvalue
      this.params.Autoencoder.muConvKernel = ones(3,3);
      this.params.Autoencoder.muFeatureKernel = [1, 0.5; 0.5 1]; % if empty then no combination in feature dimension, otherwise it should be a matrix of size (hiddenSize,hiddenSize)
      this.params.Autoencoder.muEpsilon = 1e-2;

      this.params.Autoencoder.maskingNoiseFraction = 0;
      this.params.Autoencoder.maskingNoiseFractionHidden = 0;
      this.params.Autoencoder.maskingNoiseFractionHiddenRescale = false;

      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      this.params.Autoencoder.topoNghFcn = @(x) gausswin(x,3);
      this.params.Autoencoder.topoNgh = [3 3];
      this.params.Autoencoder.topoGridDims = [4 4];
      this.params.Autoencoder.topoPeriodicBoundary = [true true];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      
      [diff,numGrad,grad,diffs] = this.checkGradient(patches);
      
    end
    
    function [diff,numGrad,grad] = unittestReLog()
      
      this = Autoencoder();
      this.params.Autoencoder.patchDimForward = 5;
      this.params.Autoencoder.patchDimBackward = 3;
      this.params.Autoencoder.hiddenSize = 16;
      this.params.Autoencoder.useSoftmax = false;
      this.params.Autoencoder.useNormMean = false;
      this.params.Autoencoder.useRectifiedLinear = false;
      this.params.Autoencoder.useRectifiedLog = true;
      this.params.Autoencoder.useNoActFcn = false;
      this.params.Autoencoder.sparsityParam = 0.035;
      this.params.Autoencoder.alpha = 1;
      this.params.Autoencoder.lambda = 0.1;
      this.params.Autoencoder.lambdaBackScale = 0;
      this.params.Autoencoder.lambdaL1 = 0;
      this.params.Autoencoder.fixL2WeightBack = true;
      this.params.Autoencoder.fixL2WeightBackPoolPerHiddenNeuron = false;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0.1;
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 2;
      this.params.Autoencoder.beta = 0;
      this.params.Autoencoder.maskingNoiseFraction = 0.25;
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      
      this.params.Autoencoder.topoNghFcn = [];
      this.params.Autoencoder.topoNgh = [1];
      this.params.Autoencoder.topoGridDims = [16];
      this.params.Autoencoder.topoPeriodicBoundary = [false];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      this.params.Autoencoder.gamma = 1;
      this.params.Autoencoder.gammaPvalue = 0.1;
      
      this.params.Autoencoder.useLinearAutoencError = false;
      patches = randn([20 20 3 3]);
      [diff,numGrad,grad] = this.checkGradient(patches);
      
    end
    
    function numgrad = computeNumericalGradient(J, theta)
      numgrad = zeros(size(theta));
      epsilon=1e-4;

      tmp = zeros(size(theta));
      
      globStream = Autoencoder.getGlobRng();
      globStreamInitState = globStream.State;

      for i=1:length(theta)
        if mod(i,100)==0
          disp(['i=' num2str(i) ' of ' num2str(length(theta))]);
        end
        tmp(max(i-1,1))=0;
        tmp(i) = 1;
        
        %Reset the random number generator always to the same value
        globStream.State = globStreamInitState;
        fvalPlus = feval(J,theta+epsilon*tmp);
        globStream.State = globStreamInitState;
        fvalMinus = feval(J,theta-epsilon*tmp);
        
        numgrad(i) = (fvalPlus - fvalMinus) / (2*epsilon);
      end
    end
    
    function rng_stream = getGlobRng()
        if sum(strcmp('getGlobalStream',methods('RandStream')))
          rng_stream = RandStream.getGlobalStream();
        else
          rng_stream = RandStream.getDefaultStream();
        end
    end
    
  end
  
end

