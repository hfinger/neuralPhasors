classdef RBM < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
    tempsavepath
    nghMat
    randstreamState
    randstreamStateBeforeReload
    randstreamStateBeforeValidation
    
    nY
    nX
    nF1
    nRfY
    nRfX
    nF2
    dimW
    
    eHid0
    eHid
    eVis
    
    W
    bforw
    bback
    dWLastGrad
    dbforwLastGrad
    dbbackLastGrad
    
    lRate
    wPenalty
    
    log = struct();
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = RBM(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.RBM.inActFolder = 'whitenedData'; %relative to the workpath
      this.params.RBM.inActFilenames = 'act.*.mat';
      this.params.RBM.outWeightsFolder = 'RBMWeights'; %relative to the workpath
      this.params.RBM.continue = false;
      this.params.RBM.inputType = 'gaussian';
      
      this.params.RBM.inSamplesDims = [64 64 3]; % [x,y,#features]
      this.params.RBM.inSamplesBorderBuffer = 0;
      
      this.params.RBM.patchDim = 8; % after subsampling of input
      this.params.RBM.hiddenSize = 100;
      this.params.RBM.inputSubsampling = 2;
      
      this.params.RBM.initWMax = 1;
      this.params.RBM.initWGaussian = false;
      this.params.RBM.initWScaleDistFcn = [];
      
      this.params.RBM.saveinterval = 10;
      this.params.RBM.saveintervalmat = 10;
      this.params.RBM.saveintervalpng = 10;
     
      this.params.RBM.topo = 0;
      this.params.RBM.topoNghFcn = [];
      this.params.RBM.topoNgh = [3 3];
      this.params.RBM.topoGridDims = [20 20];
      this.params.RBM.topoPeriodicBoundary = [true true];
      this.params.RBM.topoEpsilon = 1e-2;
      
      this.params.RBM.numberOfPatchReloads = 10;
      this.params.RBM.numberOfImagesPerPatchReload = 10;
      this.params.RBM.numberOfPatchesPerReload = 100;
      this.params.RBM.numberOfIterationsPerReload = 10;
      
      this.params.RBM.validationSetsize = 0;
      this.params.RBM.validationInterval = 10;
      
      this.params.RBM.lRate = 0.1;			% LEARNING RATE
      this.params.RBM.nGibbs = 1;				% # OF GIBBS SAMPLES (CONTRASTIVE DIVERGENCE)
      
      this.params.RBM.sparsity = 0.02;		% TARGET HIDDEN UNIT SPARSITY
      this.params.RBM.sparseGain = 1;			% GAIN ON THE LEARNING RATE FOR SPARSITY CONSTRAINTS
      this.params.RBM.momentum = 0.9;			% (DEFAULT) GRADIENT MOMENTUM FOR WEIGHT UPDATES
      this.params.RBM.wPenalty = .05;			% L2 WEIGHT PENALTY
      this.params.RBM.beginAnneal = Inf;		% BEGIN SIMULUATED ANNEALING AFTER THIS # OF EPOCHS
      this.params.RBM.beginWeightDecay = 1;	% BEGIN WEIGHT DECAY AFTER THIS # OF EPOCHS
      this.params.RBM.restrictConvPosPhase = true;
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      params = this.params.RBM;
      
      if this.numJobs > 1
        this.tempsavepath = fullfile(this.temppath,num2str(this.currJobid));
        mkdir(this.tempsavepath);
      else
        this.tempsavepath = this.temppath;
      end
      
      inputfolder = fullfile(this.workpath, params.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder, params.inActFilenames );
      
      this.nF1 = params.inSamplesDims(3)*params.inputSubsampling.^2;
      this.nY = params.inSamplesDims(1);
      this.nX = params.inSamplesDims(2);
      this.nF2 = params.hiddenSize;
      this.nRfY = params.patchDim;
      this.nRfX = params.patchDim;
      this.dimW = [params.patchDim params.patchDim this.nF1 this.nF2];
      
      if params.continue
        tmp = load(fullfile(this.tempsavepath,'currIter.mat'));
        this.W = tmp.W;
        this.bforw = tmp.bforw;
        this.bback = tmp.bback;
        this.dWLastGrad = zeros(size(this.W));
        this.dbforwLastGrad = zeros(size(this.bforw));
        this.dbbackLastGrad = zeros(size(this.bback));
        clear tmp;
      else
        if params.initWGaussian
          this.W = randn(this.dimW) * params.initWMax;
        else
          this.W = ( 2*rand(this.dimW)-1 ) * params.initWMax;
        end
        if ~isempty(params.initWScaleDistFcn)
          this.W = RBM.scaleWeightWithDist(this.W,params.initWScaleDistFcn);
        end
        this.bforw = zeros(1,1,this.nF2);
        this.bback = zeros(1,1,this.nF1);
        this.dWLastGrad = zeros(size(this.W));
        this.dbforwLastGrad = zeros(size(this.bforw));
        this.dbbackLastGrad = zeros(size(this.bback));
      end
      
      if params.topo~=0
        if isempty(params.topoNghFcn)
          this.nghMat = linearDecoderNghMat( params.topoNgh,params.topoGridDims,params.topoPeriodicBoundary);
        else
          this.nghMat = linearDecoderNghWinMat( params.topoNgh,params.topoGridDims,params.topoPeriodicBoundary, params.topoNghFcn);
        end
      end
      
      validRandStream = RandStream('mt19937ar','Seed',randi(2^32-1));
      this.randstreamStateBeforeValidation = validRandStream.State;
      
      if params.validationSetsize~=0
        if ~isempty(params.numberOfImagesPerPatchReload)
          validationSetImageIds=randi(length(pathlist),1,params.numberOfImagesPerPatchReload);
          pathlistValid = pathlist(validationSetImageIds);
          filelistValid = filelist(validationSetImageIds);
          pathlist = pathlist(setdiff(1:length(pathlist),validationSetImageIds));
          filelist = filelist(setdiff(1:length(filelist),validationSetImageIds));
          validPatches = this.preparePatches(pathlistValid,filelistValid,params.validationSetsize);
        else
          validPatches = this.preparePatches(pathlist,filelist,params.validationSetsize);
        end
      else
        validPatches = [];
      end
      
      
      for reloadCount = 1:params.numberOfPatchReloads
        
        % LOAD NEW TRAINING PATCHES
        if ~isempty(params.numberOfImagesPerPatchReload)
          imgFileIds = randi(length(pathlist),1,params.numberOfImagesPerPatchReload);
          trainPatches = this.preparePatches(pathlist(imgFileIds),filelist(imgFileIds),params.numberOfPatchesPerReload);
        else
          trainPatches = this.preparePatches(pathlist,filelist,params.numberOfPatchesPerReload);
        end
        trainPatches = trainPatches(:,:,:,randperm(size(trainPatches,4)));
        
        
        % BEGIN SIMULATED ANNEALING?
        if reloadCount < params.beginAnneal
          this.lRate = params.lRate;
        else
          this.lRate = max( params.lRate*params.beginAnneal/reloadCount, 1e-10);
        end

        % BEGIN WEIGHT DECAY?
        if reloadCount < params.beginWeightDecay
          this.wPenalty = 0;
        else
          this.wPenalty = params.wPenalty;
        end
      
        % loop several times over the currently loaded minibatch:
        for iterCount = 1:params.numberOfIterationsPerReload
          [this, minibatchlog] = this.trainBatch(trainPatches);
          fprintf('reloadCount %d --> iterCount %d --> sumErr: %f --> timeTaken: %f\n', reloadCount, iterCount,  minibatchlog.sumErr, minibatchlog.timeTaken);
          
          % log first iter:
          if iterCount==1
            this.log.trainFirstIters(reloadCount) = minibatchlog;
          end
        end
        
        % log last iter:
        this.log.trainLastIters(reloadCount) = minibatchlog;
        
        % save this class object:
        if mod(reloadCount,this.params.RBM.saveinterval)==0
          save(fullfile(this.temppath,['currIterJob' num2str(this.currJobid) '.mat']),'this');
        end
        
        % save forw conn:
        if mod(reloadCount,this.params.RBM.saveintervalpng)==0 || mod(reloadCount,this.params.RBM.saveintervalmat)==0
          
          forwConn = this.create6DConn(this.W,this.bforw);
          
          % save as conn mat:
          if mod(reloadCount,this.params.RBM.saveintervalmat)==0
            save(fullfile(this.tempsavepath,['forwConnIter' num2str(reloadCount) '.mat']),'-struct','forwConn');
          end
          
          % save as png file:
          if mod(reloadCount,this.params.RBM.saveintervalpng)==0
            totalNumFeaturesForw = size(forwConn.W,4)*size(forwConn.W,5)*size(forwConn.W,6);
            plotSizeX = round(sqrt(totalNumFeaturesForw));
            while mod(totalNumFeaturesForw/plotSizeX,1)~=0 && plotSizeX>1
              plotSizeX = plotSizeX - 1;
            end
            pngpath = fullfile(this.tempsavepath,['forwConnIter' num2str(reloadCount) '.png']);
            if size(forwConn.W,3)==1
              img =plotFeatures( reshape(forwConn.W,[size(forwConn.W,1) size(forwConn.W,2) plotSizeX totalNumFeaturesForw/plotSizeX]), [], 'gray', 2, [], [], [], false, true );
              imwrite(img,pngpath,'png');
            elseif size(forwConn.W,3)==3
              plotColorFeatures( reshape(forwConn.W,[size(forwConn.W,1) size(forwConn.W,2) size(forwConn.W,3) plotSizeX totalNumFeaturesForw/plotSizeX]), true, pngpath, true );
            end
            copyfile(pngpath,fullfile(this.temppath,['forwConnLastIterJob' num2str(this.currJobid) '.png']),'f');
          end
        end
        
      end
      
      %% save results:
      savepath = fullfile(this.workpath,params.outWeightsFolder);
      if this.numJobs > 1
        savepath = fullfile(savepath,num2str(this.currJobid));
      end
      mkdir(savepath);
      
      forwConn = this.create6DConn(this.W,this.bforw);
      backConn = this.create6DConn(permute(this.W(end:-1:1,end:-1:1,:,:),[1 2 4 3]),this.bback);
      save(fullfile(savepath,'forwConn.mat'),'-struct','forwConn')
      save(fullfile(savepath,'backConn.mat'),'-struct','backConn')
      
      disp('Finished saving connections');
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    
    function [this, log] = trainBatch(this, data)
      params = this.params.RBM;
      numSamples = size(data,4);
      log.sumErr = 0;
      sumEHid0 = 0;
      tic;
      for sampleCount = 1:numSamples
        
        % compute hidden layer expectation given the visible state:
        this = this.hidGivVis(data(:,:,:,sampleCount));
        % save it to use it later in CD learning:
        this.eHid0 = this.eHid; 
        
        % Use binary state for first hidden layer sampling (cf. chapter
        % 3.4 in "A practical guide to training restricted Boltzmann
        % machines" of G Hinton). Therefore only sample here:
        this.eHid = RBM.drawBernoulli(this.eHid); 
        
        sumEHid0 = sumEHid0 + mean(this.eHid(:));
        
        for iC = 1:params.nGibbs
          % reconstruct visible state
          this = this.visGivHid(this.eHid);
          
          % FINISH CD[n]
          this = this.hidGivVis(this.eVis);
          
          % Sample hidden, but not in last step according to Hinton Guide:
          if iC<params.nGibbs
            this.eHid = RBM.drawBernoulli(this.eHid); 
          end
        end
        
        % Contrastive Divergence Learning Rule:
        delPerSide = params.nGibbs*(params.patchDim-1);
        if params.restrictConvPosPhase
          tmpHid0 = this.eHid0(delPerSide+1:end-delPerSide,delPerSide+1:end-delPerSide,:);
          tmpData = data(delPerSide+1:end-delPerSide,delPerSide+1:end-delPerSide,:,sampleCount);
          dW = conv3dInv(tmpHid0 , tmpData);
          dW = dW - conv3dInv(this.eHid, this.eVis);
        else
          dW = conv3dInv( this.eHid0, data(:,:,:,sampleCount)) - conv3dInv(this.eHid, this.eVis);
        end
        dbback = sum(sum(data(delPerSide+1:end-delPerSide,delPerSide+1:end-delPerSide,:,sampleCount) - this.eVis,1),2);
        dbforw = sum(sum(this.eHid0(delPerSide+1:end-delPerSide,delPerSide+1:end-delPerSide,:) - this.eHid,1),2);
        
        % Add weight decay:
        if this.wPenalty~=0
          dW = dW - this.wPenalty*this.W;
        end
        
        % INDUCE HIDDEN UNIT SPARSITY (Eq 5.5; Lee, 2010)
        if params.sparseGain~=0
          dbforw = dbforw + params.sparseGain*( params.sparsity - mean(mean(this.eHid0)) );
        end
        
        % Update Model using Momentum
        [this.W, this.dWLastGrad] = this.updateParams(this.W, dW, this.dWLastGrad);
        [this.bforw, this.dbforwLastGrad] = this.updateParams(this.bforw, dbforw, this.dbforwLastGrad);
        [this.bback, this.dbbackLastGrad] = this.updateParams(this.bback, dbback, this.dbbackLastGrad);
        
        % Calc reconstruction error:
        batchErr = sum(sum(sum( ( data(delPerSide+1:end-delPerSide,delPerSide+1:end-delPerSide,:,sampleCount) - this.eVis ).^2 )));
        log.sumErr = log.sumErr + batchErr;
        
      end
      log.timeTaken = toc;
      log.meanEHid0 = sumEHid0/numSamples;
      
    end
    
    function [statevar,grad] = updateParams(this,statevar,grad,gradPrev)
      grad = this.params.RBM.momentum * gradPrev + (1-this.params.RBM.momentum)*grad;
      statevar = statevar + this.lRate * grad;
    end
    
    function this = hidGivVis(this,vis)
      this.eHid = RBM.sigmoid(bsxfun(@plus,conv3d(this.W,vis),this.bforw));
    end
    
    function this = visGivHid(this,hid)
      this.eVis = bsxfun(@plus,conv3d(this.W,hid,true),this.bback);
      if strcmp(this.params.RBM.inputType,'binary');
        this.eVis = RBM.sigmoid(this.eVis);
      end
    end
    
    function [patches] = preparePatches(this,pathlist,filelist,numSamples)
      patches = loadRandomPatches( ...
        this.params.RBM.inSamplesDims(1), ...
        this.params.RBM.inSamplesDims(2),...
        this.params.RBM.inSamplesDims(3), ...
        numSamples, ...
        pathlist, ...
        filelist, ...
        this.params.RBM.inSamplesBorderBuffer );
      
      
      if this.params.RBM.inputSubsampling > 1
        % now collapse the subsampled dimensions into the feature dim:
        patches = reshape(patches, [...
          this.params.RBM.inputSubsampling ...
          size(patches,1)/this.params.RBM.inputSubsampling ...
          this.params.RBM.inputSubsampling ...
          size(patches,2)/this.params.RBM.inputSubsampling ...
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
    
    function [forwConn] = create6DConn(this, W, b)
      
      forwConn.W = reshape(W,[size(W,1) size(W,2) size(W,3) 1 1 size(W,4)]);
      forwConn.b = b;
      %W has dims [dxin dyin fin dxout dyout fout]
      %b has dims [dxout dyout fout]
      
      if this.params.RBM.inputSubsampling > 1
        
        % inputSubsampling is so far all collapsed in fin
        % so here it is extracted and collapsed into dxin and dyin
        forwConn.W = reshape(forwConn.W,[...
          size(forwConn.W,1)...
          size(forwConn.W,2)...
          this.params.RBM.inputSubsampling...
          this.params.RBM.inputSubsampling...
          size(forwConn.W,3)/this.params.RBM.inputSubsampling^2 ...
          size(forwConn.W,4) ...
          size(forwConn.W,5) ...
          size(forwConn.W,6)]);
        forwConn.W = permute(forwConn.W,[3 1 4 2 5 6 7 8]);
        forwConn.W = reshape(forwConn.W,[...
          size(forwConn.W,1)*size(forwConn.W,2)...
          size(forwConn.W,3)*size(forwConn.W,4)...
          size(forwConn.W,5)...
          size(forwConn.W,6)...
          size(forwConn.W,7)...
          size(forwConn.W,8)]);
      end
      
      % save forward connections:
      forwConn.inputSubsampling = this.params.RBM.inputSubsampling;
      forwConn.shiftOutputdims = false; %at the moment no tiled conv implemented
      forwConn.actFcn = [];
    end
    
    function plotProgress(this)
      cmap = colormap('lines');
      
      %% Gather statistics of the following variables:
      varNames = {'eHid0','eHid','eVis', 'W', 'bforw', 'bback', 'dWLastGrad', 'dbforwLastGrad', 'dbbackLastGrad'};
      minVars = ones(this.numJobs,length(varNames))*Inf;
      maxVars = -ones(this.numJobs,length(varNames))*Inf;
      meanVar = zeros(this.numJobs,length(varNames));
      
      %% Job names:
      jobNames = cell(1,this.numJobs);
      for jobid=1:this.numJobs
        jobNames{jobid} = '';
        for paramid=1:size(this.paramComb,1)
          jobNames{jobid} = [jobNames{jobid} ' ' this.variableParams{paramid}{2} '=' num2str(this.paramComb{paramid,jobid})];
        end
      end
      
      %% plot specific log variable:
      for jobid=1:this.numJobs
        currIter = load(fullfile(this.temppath,['currIterJob' num2str(jobid) '.mat']),'this');
        
        % gather statistics of min, mean and max for later plots:
        for varid=1:length(varNames)
          minVars(jobid,varid)=min(currIter.this.(varNames{varid})(:));
          maxVars(jobid,varid)=max(currIter.this.(varNames{varid})(:));
          meanVar(jobid,varid)=mean(currIter.this.(varNames{varid})(:));
        end
        
        if jobid==1
          fnames = fieldnames(currIter.this.log.trainFirstIters);
          [Selection,ok] = listdlg('PromptString','Show logs of the following variables:','ListString',fnames);
          if ok
            figure(1);
          end
        end
        
        if ok
          plot([currIter.this.log.trainFirstIters.(fnames{Selection})],'Color',cmap(jobid,:));
          hold on;
        end
      end
      if ok
        hold off;
        legend(jobNames(:))
      end
      
      %% display min, max and mean of variables:
      [varid,ok] = listdlg('PromptString','Show min,max and mean of the following variables:','ListString',varNames);
      if ok
        fprintf('Statistics of %s\n', varNames{varid});
        for jobid=1:this.numJobs
          fprintf('job: %s --> min: %f --> mean: %f --> max: %f\n', jobNames{jobid}, minVars(jobid,varid), meanVar(jobid,varid), maxVars(jobid,varid));
        end
      end
      
      %% plot hist:
      [varid,ok] = listdlg('PromptString','Plot Histograms of the following variables:','ListString',varNames);
      if ok
        if Selection==1
          figure(2);
          edges = linspace(min(minVars(:,varid)),max(maxVars(:,varid)),100); 
          for jobid=1:this.numJobs
            currIter = load(fullfile(this.temppath,['currIterJob' num2str(jobid) '.mat']),'this');
            n = histc(currIter.this.(varNames{varid})(:),edges);
            plot(edges,n,'Color',cmap(jobid,:));
            hold on;
          end
          hold off;
          legend(jobNames(:))
        end
      end
    end
    
  end
  
  methods (Static)
    
    function W = scaleWeightWithDist(W,distInitFcn)
      distMat = sqrt(bsxfun(@plus,((1:size(W,1))-(size(W,1)+1)/2).^2', ((1:size(W,2))-(size(W,2)+1)/2).^2));
      distMat = feval(distInitFcn,distMat);
      distMat = distMat / mean(distMat(:));
      W=bsxfun(@times,W,distMat);
    end
    
    function sigm = sigmoid(x)
      sigm = 1 ./ (1 + exp(-x));
    end
    
    function p = drawBernoulli(p)
      % p = drawNormal(mu);
      %--------------------------------------------------------------------------
      % Draw samples from a multivariate normal  with mean <mu> and identity
      % covariance.
      %--------------------------------------------------------------------------
      p = double(rand(size(p)) < p);
    end
    
  end
  
  
end

