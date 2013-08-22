classdef Autoencoder < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
    tempsavepath
    forwardRfSize
    backwardRfSize
    nghMat
    randstreamState
    randstreamStateBeforeReload
    randstreamStateBeforeValidation
    
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
      this.params.Autoencoder.outWeightsFolder = 'AEWeights'; %relative to the workpath
      this.params.Autoencoder.continue = false;
      this.params.Autoencoder.debugOutput = false;
            
      this.params.Autoencoder.continueBatch = false; % !!!!!!! not working correctly! TODO: reorder dimensions of W!
      this.params.Autoencoder.continueBatchInWeightsFolder = 'AEWeights';
      this.params.Autoencoder.continueBatchInBackConnFilenames = 'backConnIter50.mat';
      this.params.Autoencoder.continueBatchInForwConnFilenames = 'forwConnIter50.mat';
      
      this.params.Autoencoder.loadPatchesAndTheta = [];
      
      this.params.Autoencoder.inSamplesDims = [32 32 3 500]; % [x,y,#features,#samples]
      this.params.Autoencoder.inSamplesBorderBuffer = 0;
      
      this.params.Autoencoder.patchDimForward = 8;
      this.params.Autoencoder.patchDimBackward = 8;
      this.params.Autoencoder.hiddenSize = 400;
      this.params.Autoencoder.inputSubsampling = 1;
      
%       this.params.Autoencoder.forwDims = [12 12 3 1 1 1 1 100]; %xoutTile,youtTile,dxoutShift,dyoutShift,dxout,dyout,fout
%       this.params.Autoencoder.forwStride = 1;
%       this.params.Autoencoder.forwShift = false;
      this.params.Autoencoder.forwInitMaxWeight = 1;
      this.params.Autoencoder.forwInitScaleDistFcn = [];
      
%       this.params.Autoencoder.backDims = [12 12 100 1 1 1 1 3];
%       this.params.Autoencoder.backStride = 1;
%       this.params.Autoencoder.backShift = false;
      this.params.Autoencoder.backInitMaxWeight = 1;
      this.params.Autoencoder.backInitScaleDistFcn = [];
      
%       this.params.Autoencoder.symWeights = false;
      
      this.params.Autoencoder.useRectifiedLinear = false; %instead of sigmoid
      this.params.Autoencoder.useNoActFcn = false; %instead of sigmoid
      this.params.Autoencoder.useSoftmax = false;
%       this.params.Autoencoder.conv3tiledMethod = 'prod'; %choose between 'prod','tprod','im2col'
      
      this.params.Autoencoder.saveinterval = 10;
      this.params.Autoencoder.saveintervalpng = 10;
      this.params.Autoencoder.savepng = false;
      this.params.Autoencoder.sparsityParam = 0.035; % desired average activation of the hidden units.
      
      this.params.Autoencoder.alpha = 1;        % weight of Autoencoder
      this.params.Autoencoder.lambda = 3e-3;    % weight of L2-decay parameter
      this.params.Autoencoder.lambdaL1 = 0 ;    % weight of L1-decay parameter
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0 ; %if not zero then smooth it
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 0 ;
      this.params.Autoencoder.beta = 5;         % weight of sparsity penalty term
      this.params.Autoencoder.gamma = 0;        % weight of topographic sparsity penalty term
      this.params.Autoencoder.maskingNoiseFraction = 0;        % fraction of inputs to set to 0
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      
      this.params.Autoencoder.topoNghFcn = [];
      this.params.Autoencoder.topoNgh = [3 3];
      this.params.Autoencoder.topoGridDims = [20 20];
      this.params.Autoencoder.topoPeriodicBoundary = [true true];
      this.params.Autoencoder.topoEpsilon = 1e-2;
    
      this.params.Autoencoder.numberOfPatchReloads = 1;
      this.params.Autoencoder.numberOfImagesPerPatchReload = [];
      
      this.params.Autoencoder.batchsize = []; %[] means full batch
      this.params.Autoencoder.fixedBatches = false;
      
      this.params.Autoencoder.validationSetsize = 0;
      this.params.Autoencoder.validationInterval = 10;
      this.params.Autoencoder.validationSetIds = [];
      this.params.Autoencoder.validationSetImageIds = []; % will only be used if also numberOfImagesPerPatchReload is set
      
      this.params.Autoencoder.resetRandstreamEachReload = false;
      this.params.Autoencoder.resetRandstreamEachIter = false;
      this.params.Autoencoder.resetRandstreamEachEval = false;
      
      this.params.Autoencoder.useMinFuncGrad = true;
      
      this.params.minFuncGrad.Method = 'rmsprop';
      this.params.minFuncGrad.maxIter = 500;
      this.params.minFuncGrad.learnrate = 1e-5;
      this.params.minFuncGrad.EMAconst = 0.1;
      this.params.minFuncGrad.momentum = 0.9; % = 1 - 2/(lambda*2.8854+1) for halflife lambda
      this.params.minFuncGrad.display = 'on';
      
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
      
      if this.numJobs > 1
        this.tempsavepath = fullfile(this.temppath,num2str(this.currJobid));
        mkdir(this.tempsavepath);
      else
        this.tempsavepath = this.temppath;
      end
      
      inputfolder = fullfile(this.workpath, this.params.Autoencoder.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder, this.params.Autoencoder.inActFilenames );
      
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
      else
        validPatches = [];
      end
      
      this.forwardRfSize = [this.params.Autoencoder.patchDimForward this.params.Autoencoder.patchDimForward this.params.Autoencoder.inSamplesDims(3)*this.params.Autoencoder.inputSubsampling.^2];
      this.backwardRfSize = [this.params.Autoencoder.patchDimBackward this.params.Autoencoder.patchDimBackward this.params.Autoencoder.hiddenSize];
      dimW1 = [this.forwardRfSize 1 1 this.backwardRfSize(3)];
      dimW2 = [this.backwardRfSize 1 1 this.forwardRfSize(3)];
      
      if this.params.Autoencoder.continueBatch
        inWpath = fullfile(this.workpath,this.params.Autoencoder.continueBatchInWeightsFolder);
        if this.numJobs > 1
          inWpath = fullfile(inWpath,num2str(this.currJobid));
        end
        tmp = load(fullfile(inWpath,this.params.Autoencoder.continueBatchInForwConnFilenames));
        % !!!!!!! not working correctly! TODO: reorder dimensions of W!
        W1 = reshape(tmp.W,dimW1);
        b1 = tmp.b;
        tmp = load(fullfile(inWpath,this.params.Autoencoder.continueBatchInBackConnFilenames));
        % !!!!!!! not working correctly! TODO: reorder dimensions of W!
        W2 = reshape(tmp.W,dimW2);
        b2 = tmp.b;
        clear tmp;
      elseif this.params.Autoencoder.continue
        tmp = load(fullfile(this.tempsavepath,'currIter.mat'));
        W1 = tmp.W1;
        W2 = tmp.W2;
        b1 = tmp.b1;
        b2 = tmp.b2;
        clear tmp;
      else
        W1 = ( 2*rand(dimW1)-1 ) * this.params.Autoencoder.forwInitMaxWeight;
        W2 = ( 2*rand(dimW2)-1 ) * this.params.Autoencoder.backInitMaxWeight;
        if ~isempty(this.params.Autoencoder.forwInitScaleDistFcn)
          W1 = Autoencoder.scaleWeightWithDist(W1,this.params.Autoencoder.forwInitScaleDistFcn);
        end
        if ~isempty(this.params.Autoencoder.backInitScaleDistFcn)
          W2 = Autoencoder.scaleWeightWithDist(W2,this.params.Autoencoder.backInitScaleDistFcn);
        end
        b1 = zeros(dimW1(4:end));
        b2 = zeros(dimW2(4:end));
      end
      theta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
      
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
      
      validRandStream = RandStream('mt19937ar','Seed',randi(2^32-1));
      this.randstreamStateBeforeValidation = validRandStream.State;
      
      globStream = RandStream.getDefaultStream;
      if this.params.Autoencoder.resetRandstreamEachReload
        this.randstreamStateBeforeReload = globStream.State;
      end
      
      for reloadId = 0:this.params.Autoencoder.numberOfPatchReloads-1
        numTrainPatches = this.params.Autoencoder.inSamplesDims(4) - this.params.Autoencoder.validationSetsize;
        if ~isempty(this.params.Autoencoder.numberOfImagesPerPatchReload)
          imgFileIds = randi(length(pathlist),1,this.params.Autoencoder.numberOfImagesPerPatchReload);
          if this.params.Autoencoder.debugOutput
            disp(['New selected imgFileIds: ' num2str(imgFileIds(1:4))])
          end
          patches = this.preparePatches(pathlist(imgFileIds),filelist(imgFileIds),numTrainPatches);
        else
          patches = this.preparePatches(pathlist,filelist,numTrainPatches);
        end
        
        if this.params.Autoencoder.resetRandstreamEachReload
          globStream.State = this.randstreamStateBeforeReload;
        end
        
        if this.params.Autoencoder.resetRandstreamEachIter || this.params.Autoencoder.resetRandstreamEachEval
          this.randstreamState = globStream.State;
        end
        
        if this.params.Autoencoder.useMinFuncGrad
          options = this.params.minFuncGrad;
          options.outputFcn = @(optTheta,iterationType,i,funEvals,f,t,gtd,g,d,optCond) this.iterFcnSave( optTheta,iterationType,i+reloadId*options.maxIter,funEvals,f,t,gtd,g,d,optCond, validPatches);
  %         save(fullfile(this.tempsavepath,'initialParams.mat'),'-v7.3','theta','this','patches','options');%
          theta = minFuncGrad( @(p) autoencoderCost(this, p, patches), theta, options);
        else
          options = this.params.minFunc;
          options.outputFcn = @(optTheta,iterationType,i,funEvals,f,t,gtd,g,d,optCond) this.iterFcnSave( optTheta,iterationType,i+reloadId*options.maxIter,funEvals,f,t,gtd,g,d,optCond, validPatches);
  %         save(fullfile(this.tempsavepath,'initialParams.mat'),'-v7.3','theta','this','patches','options');
          theta = minFunc( @(p) autoencoderCost(this, p, patches), theta, options);
        end
        
      end
      
      %% save results:
      savepath = fullfile(this.workpath,this.params.Autoencoder.outWeightsFolder);
      if this.numJobs > 1
        savepath = fullfile(savepath,num2str(this.currJobid));
      end
      mkdir(savepath);
      
      [forwConn, backConn] = this.createConn(theta); %#ok<ASGLU,NASGU>
      save(fullfile(savepath,'forwConn.mat'),'-struct','forwConn')
      save(fullfile(savepath,'backConn.mat'),'-struct','backConn')
      
      disp('Finished saving connections');
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function [patches] = preparePatches(this,pathlist,filelist,numSamples)
      patches = loadRandomPatches( ...
        this.params.Autoencoder.inSamplesDims(1), ...
        this.params.Autoencoder.inSamplesDims(2),...
        this.params.Autoencoder.inSamplesDims(3), ...
        numSamples, ...
        pathlist, ...
        filelist, ...
        this.params.Autoencoder.inSamplesBorderBuffer );
      
      
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
      elseif this.params.Autoencoder.useNoActFcn
        forwConn.actFcn = [];
      else
        forwConn.actFcn = @(x) 1 ./ (1 + exp(-x));
      end
      
      if this.params.Autoencoder.useSoftmax
        forwConn.actFcn2 = @(x) feval(@(expx) bsxfun(@rdivide,expx,sum(expx,3)), exp(x));
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
      
      if mod(i,this.params.Autoencoder.saveintervalpng)==0 || mod(i,this.params.Autoencoder.saveinterval)==0
        [forwConn, backConn] = this.createConn(optTheta); %#ok<NASGU>
      end
      
      if mod(i,this.params.Autoencoder.saveinterval)==0
        save(fullfile(this.tempsavepath,['forwConnIter' num2str(i) '.mat']),'-struct','forwConn');
        save(fullfile(this.tempsavepath,['backConnIter' num2str(i) '.mat']),'-struct','backConn');
        save(fullfile(this.tempsavepath,['trainInfoIter' num2str(i) '.mat']),'iterationType','i','funEvals','f','t','gtd','d','optCond');
      end
      if mod(i,this.params.Autoencoder.saveintervalpng)==0
        if this.params.Autoencoder.savepng
          totalNumFeaturesForw = size(forwConn.W,4)*size(forwConn.W,5)*size(forwConn.W,6);
          pngpath = fullfile(this.tempsavepath,['forwConnIter' num2str(i) '.png']);
          plotColorFeatures(  reshape(forwConn.W,[size(forwConn.W,1) size(forwConn.W,2) size(forwConn.W,3) sqrt(totalNumFeaturesForw) sqrt(totalNumFeaturesForw)]), true, pngpath, true );
          copyfile(pngpath,fullfile(this.temppath,['forwConnLastIterJob' num2str(this.currJobid) '.png']),'f');
        end
      end
      
      %% do validation set:
      if this.params.Autoencoder.validationSetsize~=0
        if mod(i,this.params.Autoencoder.validationInterval)==0
          
          % set both randstream states such that it will be used during validation:
          globStream = RandStream.getDefaultStream;
          randstreamStateTmp1 = globStream.State;
          globStream.State = this.randstreamStateBeforeValidation;
          randstreamStateTmp2 = this.randstreamState;
          this.randstreamState = this.randstreamStateBeforeValidation;
          
          if ~isempty(this.params.Autoencoder.batchsize)
            tmp = minibatchIds;
            minibatchIds = 1:this.params.Autoencoder.validationSetsize;
          end
          
          [cost,~,optout] = this.autoencoderCost( optTheta, validPatches);
          disp(['iter ' num2str(i) ' validation cost: ' num2str(cost)])
          iterationType = struct();
          iterationType.optout = optout;
          iterationType.cost = cost;
          iterationType.iteration = i; %#ok<STRNU>
          save(fullfile(this.tempsavepath,['validInfoIter' num2str(i) '.mat']),'iterationType','i','funEvals','f','t','gtd','d','optCond');          
          if ~isempty(this.params.Autoencoder.batchsize)
            minibatchIds = tmp;
          end
          
          % reset both randstream states back to the state before validation:
          globStream.State = randstreamStateTmp1;
          this.randstreamState = randstreamStateTmp2;
          
        end
      end
      
      this.selectMinibatch();
      
      if this.params.Autoencoder.resetRandstreamEachIter
        globStream = RandStream.getDefaultStream;
        globStream.State = this.randstreamState;
      end
      
      stop=false;
    end
    
    
    function [W1, W2, b1, b2] = decodeTheta( this, theta )
      
      dimW1 = [this.forwardRfSize(1), this.forwardRfSize(2), this.forwardRfSize(3), 1, 1, this.backwardRfSize(3)];
      dimW2 = [this.backwardRfSize(1), this.backwardRfSize(2), this.backwardRfSize(3), 1, 1, this.forwardRfSize(3)];
      dimb1 = [1, 1, this.backwardRfSize(3)];
      dimb2 = [1, 1, this.forwardRfSize(3)];
      
      W1 = theta(1:prod(dimW1));
      W2 = theta(prod(dimW1)+1:prod(dimW1)+prod(dimW2));
      b1 = theta(prod(dimW1)+prod(dimW2)+1:prod(dimW1)+prod(dimW2)+prod(dimb1));
      b2 = theta(prod(dimW1)+prod(dimW2)+prod(dimb1)+1:prod(dimW1)+prod(dimW2)+prod(dimb1)+prod(dimb2));
      
      W1 = reshape(W1,dimW1);
      W2 = reshape(W2,dimW2);
      b1 = reshape(b1,dimb1);
      b2 = reshape(b2,dimb2);
      
    end
    
    
    %% This is the cost function:
    function [cost,grad,optout] = autoencoderCost(this, theta, data)
      global minibatchIds;
      
      if this.params.Autoencoder.resetRandstreamEachEval
        globStream = RandStream.getDefaultStream;
        globStream.State = this.randstreamState;
      end
      
      if this.params.Autoencoder.debugOutput
        tmpglobStream = RandStream.getDefaultStream;
        tmpglobStreamInitState = tmpglobStream.State;
        disp(['globstream state: ' num2str(tmpglobStreamInitState(1:4)')])
      end
      
      if ~isempty(this.params.Autoencoder.batchsize)
        data=data(:,:,:,minibatchIds);
      end
      params = this.params.Autoencoder;
      
      [W1, W2, b1, b2] = this.decodeTheta( theta );
      % Now: W has dims [dxin dyin fin dxout dyout fout]
      % Now: b has dims [dxout dyout fout]
      W1 = permute(W1,[1 2 3 6 4 5]);
      W2 = permute(W2,[1 2 3 6 4 5]);
      b1 = permute(b1,[3 1 2]);
      b2 = permute(b2,[3 1 2]);
      % Now: W has dims [dxin dyin fin fout 1 1]
      % Now: b has dims [fout 1 1]
      % Remember that inputSubsampling is all collapsed in fin

      m=size(data,4);
      
      if params.beta ~= 0
        
        if params.maskingNoiseFraction || params.gaussianNoiseSigmaInput
          globStream = RandStream.getDefaultStream;
          globStreamInitState = globStream.State;
        end
        
        avgAct2counter = 0;
        avgAct2=zeros([1 1 numel(b1)]);
        for k=1:m
          
          curData = data(:,:,:,k);

          if params.maskingNoiseFraction
            tmp=randperm(numel(curData));
            numelToSetZero = round(params.maskingNoiseFraction*numel(curData));
            curData(tmp(1:numelToSetZero)) = 0;
            realRatioSetZero = numelToSetZero/numel(curData);
            curData = curData / (1-realRatioSetZero);
          end
          if params.gaussianNoiseSigmaInput
            curData = curData + params.gaussianNoiseSigmaInput * randn(size(curData));
          end

          z2 = conv3d(W1,curData);
          z2 = bsxfun(@plus,z2,reshape(b1,[1 1 numel(b1)]));
          
          if params.useRectifiedLinear
            a2 = max(z2,0);
          elseif params.useNoActFcn
            a2 = z2;
          else
            a2 = Autoencoder.sigmoid(z2);
          end
          
%           if params.useSoftmax
%             y2 = Autoencoder.softmax(a2);
%           end
          
          avgAct2 = avgAct2 + sum(sum(a2,2),1);
          avgAct2counter = avgAct2counter + size(a2,1)*size(a2,2);
        end
        avgAct2 = avgAct2/avgAct2counter;
        
        if params.maskingNoiseFraction || params.gaussianNoiseSigmaInput
          globStream.State = globStreamInitState;
        end
      end
      
      costAutoenc = 0;
      costRegularize = (params.lambda / 2) * (sum(W1(:).^2) + sum(W2(:).^2));
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
            costRegularizeL1 = costRegularizeL1 + params.lambdaL1 * sum(distMatW1(:) .* sum(reshape(W1L1,[size(W1,1)*size(W1,2) size(W1,3)*size(W1,4)*size(W1,5)*size(W1,5)]),2));
            costRegularizeL1 = costRegularizeL1 + params.lambdaL1 * sum(distMatW2(:) .* sum(reshape(W2L1,[size(W2,1)*size(W2,2) size(W2,3)*size(W2,4)*size(W2,5)*size(W2,6)]),2));
          else
            costRegularizeL1 = params.lambdaL1 * (sum(W1L1(:)) + sum(W2L1(:)));
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
      
      if nargout>2
        sumSquaredDelta2Autoenc = 0;
        sumSquaredDelta2Sparse = 0;
        sumSquaredDelta2Toposparse = 0;
      end
      
      for k=1:m
        curData = data(:,:,:,k);
        
        if params.maskingNoiseFraction
          tmp=randperm(numel(curData));
          numelToSetZero = round(params.maskingNoiseFraction*numel(curData));
          curData(tmp(1:numelToSetZero)) = 0;
          realRatioSetZero = numelToSetZero/numel(curData);
          curData = curData / (1-realRatioSetZero);
        end
        if params.gaussianNoiseSigmaInput
          curData = curData + params.gaussianNoiseSigmaInput * randn(size(curData));
        end
        
        z2 = conv3d(W1,curData);
        z2=bsxfun(@plus,z2,reshape(b1,[1 1 numel(b1)]));
        if params.useRectifiedLinear
          a2=max(z2,0);
        elseif params.useNoActFcn
          a2 = z2;
        else
          a2 = Autoencoder.sigmoid(z2);
        end
        
        if params.useSoftmax
          y2 = Autoencoder.softmax(a2);
        end
        
        if params.gamma~=0
          a2reduced = reshape(a2,[size(a2,1)*size(a2,2) size(a2,3)]);
          Psi_sparse = sqrt(a2reduced.^2 * this.nghMat + params.topoEpsilon);
          avgPsi_sparse = avgPsi_sparse + sum(Psi_sparse,1);
        end
        
        if params.useSoftmax
          z3 = conv3d(W2,y2);
        else
          z3 = conv3d(W2,a2);
        end
        z3 = bsxfun(@plus,z3,reshape(b2,[1 1 numel(b2)]));
        a3 = z3;
        
        %now cut out the center of the input data to match layer a3
        cut1=(size(data,1)-size(a3,1))/2;
        cut2=(size(data,2)-size(a3,2))/2;
        costAutoenc = costAutoenc + sum(reshape(a3 - data(cut1+1:end-cut1,cut2+1:end-cut2,:,k),[1 numel(a3)]).^2 / 2);
        
        if nargout > 1
          delta3 = - (data(cut1+1:end-cut1,cut2+1:end-cut2,:,k) - a3) * params.alpha;
          delta3Full = zeros(size(delta3)+[(this.backwardRfSize(1)-1)*2 (this.backwardRfSize(2)-1)*2 0]);%zeros(size(data(:,:,:,k)));
          cut1Delta3Full = (this.backwardRfSize(1)-1);
          cut2Delta3Full = (this.backwardRfSize(2)-1);
          delta3Full(cut1Delta3Full+1:end-cut1Delta3Full,cut2Delta3Full+1:end-cut2Delta3Full,:) = delta3;
          
          delta2 = zeros(size(a2));
          if params.alpha~=0
            delta2 = bsxfun(@plus,delta2,conv3d(permute(flipdim(flipdim(W2,1),2),[1 2 4 3]),delta3Full));
            if nargout>2
              sumSquaredDelta2Autoenc = sumSquaredDelta2Autoenc + sum(delta2(:).^2);
            end
          end
          
          % Here is delta2=dPsi_da2 
          if params.useSoftmax
            % Here is delta2=dPsi_da2 where a2=exp(y2)/sum(exp(y2),3)
            % where y2 is the activation after sigmoid or linear threshold
            tmpDelta2a2 = delta2 .* y2;
            
            % In several lines of code:
            %tmpDelta2a2SumWithoutOne = bsxfun(@minus, sum( tmpDelta2a2, 3), tmpDelta2a2);
            %tmpDelta2smax_notSameCoord = tmpDelta2a2SumWithoutOne .* a2;
            %tmpDelta2smax_sameCoord = delta2 .* (a2.*(1-a2));
            %delta2 = tmpDelta2smax_sameCoord - tmpDelta2smax_notSameCoord;
            
            % The same in one line of code:
            delta2 = delta2 .* (y2.*(1-y2)) - bsxfun(@minus, sum( tmpDelta2a2, 3), tmpDelta2a2) .* y2;
            
            %a2=y2; %go back to variable names as if no softmax
          end
          
          if params.beta~=0
            delta2Sparse = params.beta*(- params.sparsityParam./avgAct2 + (1-params.sparsityParam)./(1-avgAct2) )  / numUnits2;
            delta2 = bsxfun(@plus,delta2,delta2Sparse);
            if nargout>2
              sumSquaredDelta2Sparse = sumSquaredDelta2Sparse + sum(delta2Sparse(:).^2);
            end
          end
          if params.gamma~=0
            deltaSparse = zeros(size(delta2,1)*size(delta2,2),size(delta2,3));
            for f2=1:size(z2,3)
              %TODO: change to im2col:
              dPsiToposparse_da2 = bsxfun(@times,a2reduced(:,f2),this.nghMat(f2,:)) ./ Psi_sparse;
              deltaSparse(:,f2) = deltaSparse(:,f2) + sum( dPsiToposparse_da2, 2 );
            end
            delta2Toposparse = params.gamma*reshape(deltaSparse,size(delta2)) / numUnits2;
            delta2 = bsxfun(@plus,delta2, delta2Toposparse);
            if nargout>2
              sumSquaredDelta2Toposparse = sumSquaredDelta2Toposparse + sum(delta2Toposparse(:).^2);
            end
          end
          
          
          
          if params.useRectifiedLinear
            delta2 = delta2 .* sign(a2);
          elseif params.useNoActFcn
            %do nothing
          else
            delta2 = delta2 .* (a2.*(1-a2));
          end
          % Here is delta2=dPsi_dz2
          
          W1grad = W1grad + Autoencoder.conv3d3(delta2,curData);
          if params.useSoftmax
            W2grad = W2grad + Autoencoder.conv3d3(delta3,y2);
          else
            W2grad = W2grad + Autoencoder.conv3d3(delta3,a2);
          end
          
          b1grad = b1grad + squeeze(sum(sum(delta2,2),1));
          b2grad = b2grad + squeeze(sum(sum(delta3,2),1));
        end
      end
      
      if params.gamma~=0
        avgPsi_sparse = avgPsi_sparse/(m*numUnits2);
        costTopoSparse = params.gamma*sum(avgPsi_sparse);
      else
        costTopoSparse = 0;
      end
      costAutoenc = params.alpha*costAutoenc / m;
      b1grad = b1grad / m;
      b2grad = b2grad / m;
      W1grad = W1grad / m;
      W2grad = W2grad / m;
      if nargout>2
        optout.sumSquaredW1gradAutoencSparse = sum(W1grad(:).^2);
        optout.sumSquaredW2gradAutoencSparse = sum(W2grad(:).^2);
      end
      
      if params.lambda~=0
        if nargout>2
          optout.sumSquaredW1gradRegL2 = params.lambda * sum(W1(:).^2);
          optout.sumSquaredW2gradRegL2 = params.lambda * sum(W2(:).^2);
        end
        W1grad = W1grad + params.lambda * W1;
        W2grad = W2grad + params.lambda * W2;
      end
      
      if params.lambdaL1~=0
        if params.smoothLambdaL1WithEpsilon>0
          W1gradRegL1 = params.lambdaL1 * W1./sqrt(params.smoothLambdaL1WithEpsilon+W1.^2);
          W2gradRegL1 = params.lambdaL1 * W2./sqrt(params.smoothLambdaL1WithEpsilon+W2.^2);
        else
          W1gradRegL1 = params.lambdaL1 * sign(W1);
          W2gradRegL1 = params.lambdaL1 * sign(W2);
        end
        if params.scaleLambdaL1WithDistExponent
          W1gradRegL1 = bsxfun(@times,distMatW1,W1gradRegL1);
          W2gradRegL1 = bsxfun(@times,distMatW2,W2gradRegL1);
        end
        W1grad = W1grad + W1gradRegL1;
        W2grad = W2grad + W2gradRegL1;
        if nargout>2
          optout.sumSquaredW1gradRegL1 = sum(W1gradRegL1(:).^2);
          optout.sumSquaredW2gradRegL1 = sum(W2gradRegL1(:).^2);
        end
      end
      
      if nargout > 1
        grad = [W1grad(:) ; W2grad(:) ; b1grad(:) ; b2grad(:)];
      end
      cost = costAutoenc + costRegularize + costRegularizeL1 + costSparseness + costTopoSparse;
      
      if nargout>2
        optout.costAutoenc = costAutoenc;
        optout.costRegularize = costRegularize;
        optout.costRegularizeL1 = costRegularizeL1;
        optout.costSparseness = costSparseness;
        optout.costTopoSparse = costTopoSparse;
        %   optout.grad.b1grad = b1grad;
        %   optout.grad.b2grad = b2grad;
        %   optout.grad.W1grad = W1grad;
        %   optout.grad.W2grad = W2grad;
        
        if params.alpha~=0
          optout.sumSquaredDelta2Autoenc =sumSquaredDelta2Autoenc;
        end
        if params.beta~=0
          optout.sumSquaredDelta2Sparse =sumSquaredDelta2Sparse;
        end
        if params.gamma~=0
          optout.sumSquaredDelta2Toposparse =sumSquaredDelta2Toposparse;
        end
      end
      
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

      if this.numJobs > 1
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
      else
%         weights=load(fullfile(this.temppath,'currIter.mat'),'W1');
%         if exist('whitening','var')
%           weights.W1 = addWhiteningWeights(weights.W1, whitening.W);
%         end
%         savepath = fullfile(this.resultpath, 'weights.png');
%         plotColorFeatures(  reshape(weights.W1,[12 12 3 10 size(weights.W1,4)/10]), true, savepath, true );
      end
      
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
    
    
    function [diff,numGrad,grad] = checkGradient(this, patches)
      
      this.forwardRfSize = [this.params.Autoencoder.patchDimForward this.params.Autoencoder.patchDimForward size(patches,3)];
      this.backwardRfSize = [this.params.Autoencoder.patchDimBackward this.params.Autoencoder.patchDimBackward this.params.Autoencoder.hiddenSize];
      
      W1 = randn([this.forwardRfSize this.backwardRfSize(3)]);
      W2 = randn([this.backwardRfSize this.forwardRfSize(3)]);
      b1 = randn(this.backwardRfSize(3),1);
      b2 = randn(this.forwardRfSize(3),1);
      theta = [W1(:) ; W2(:) ; b1(:) ; b2(:)];
      
      if this.params.Autoencoder.gamma~=0
        if isempty(this.params.Autoencoder.topoNghFcn)
          this.nghMat = linearDecoderNghMat( this.params.Autoencoder.topoNgh,this.params.Autoencoder.topoGridDims,this.params.Autoencoder.topoPeriodicBoundary);
        else
          this.nghMat = linearDecoderNghWinMat( this.params.Autoencoder.topoNgh,this.params.Autoencoder.topoGridDims,this.params.Autoencoder.topoPeriodicBoundary, this.params.Autoencoder.topoNghFcn);
        end
      end
      
      globStream = RandStream.getDefaultStream;
      globStreamInitState = globStream.State;
      [~, grad] = this.autoencoderCost(theta, patches);
      globStream.State = globStreamInitState;
      numGrad = Autoencoder.computeNumericalGradient(@(p) autoencoderCost(this, p, patches), theta);
      
      diff = norm(numGrad-grad)/norm(numGrad+grad);
      
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
      this.params.Autoencoder.lambda = 1;
      this.params.Autoencoder.lambdaL1 = 1;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0.1;
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 2;
      this.params.Autoencoder.beta = 1;
      this.params.Autoencoder.gamma = 1;
      this.params.Autoencoder.maskingNoiseFraction = 0.25;
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0.5;
      this.params.Autoencoder.topoNghFcn = @(x) gausswin(x,3);
      this.params.Autoencoder.topoNgh = [3 3];
      this.params.Autoencoder.topoGridDims = [4 4];
      this.params.Autoencoder.topoPeriodicBoundary = [true true];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      patches = randn([20 20 3 3]);
      [diff,numGrad,grad] = this.checkGradient(patches);
      disp(['normalized diff grad: ' num2str(diff)]); 
      disp([numGrad(1:5),grad(1:5)])
      
    end
    
    function unittestReLU()
      
      this = Autoencoder();
      this.params.Autoencoder.patchDimForward = 5;
      this.params.Autoencoder.patchDimBackward = 3;
      this.params.Autoencoder.hiddenSize = 16;
      this.params.Autoencoder.useSoftmax = true;
      this.params.Autoencoder.useRectifiedLinear = true;
      this.params.Autoencoder.useNoActFcn = true;
      this.params.Autoencoder.sparsityParam = 0.035;
      this.params.Autoencoder.alpha = 1;
      this.params.Autoencoder.lambda = 0.1;
      this.params.Autoencoder.lambdaL1 = 0;
      this.params.Autoencoder.smoothLambdaL1WithEpsilon = 0.1;
      this.params.Autoencoder.scaleLambdaL1WithDistExponent = 2;
      this.params.Autoencoder.beta = 0;
      this.params.Autoencoder.gamma = 1;
      this.params.Autoencoder.maskingNoiseFraction = 0.25;
      this.params.Autoencoder.gaussianNoiseSigmaInput = 0;
      this.params.Autoencoder.topoNghFcn = [];
      this.params.Autoencoder.topoNgh = [1];
      this.params.Autoencoder.topoGridDims = [16];
      this.params.Autoencoder.topoPeriodicBoundary = [false];
      this.params.Autoencoder.topoEpsilon = 1e-2;
      patches = randn([20 20 3 3]);
      [diff,numGrad,grad] = this.checkGradient(patches);
      disp(['normalized diff grad: ' num2str(diff)]); 
      disp([numGrad(1:5),grad(1:5)])
      
    end
    
    
    function numgrad = computeNumericalGradient(J, theta)
      numgrad = zeros(size(theta));
      epsilon=1e-4;

      tmp = zeros(size(theta));
      
      globStream = RandStream.getDefaultStream;
      globStreamInitState = globStream.State;

      for i=1:length(theta)
%         disp(['i=' num2str(i) ' of ' num2str(length(theta))]);
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
    
  end
  
end

