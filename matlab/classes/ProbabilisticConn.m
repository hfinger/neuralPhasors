classdef ProbabilisticConn < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ProbabilisticConn(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ProbabilisticConn.inCorrFile = 'FeatureCorrCov/patchCorr.mat'; %relative to the workpath
      this.params.ProbabilisticConn.outWeightsFolder = 'SparseIntralayerConn'; %relative to the workpath
      this.params.ProbabilisticConn.numExc = 1000;
      this.params.ProbabilisticConn.numInh = 200;
      this.params.ProbabilisticConn.exclSelfConn = true;
      this.params.ProbabilisticConn.useDiscretesample = false; %% better set to true ?
      this.params.ProbabilisticConn.exclDoubleConns = false;
      this.params.ProbabilisticConn.probRadiusFcn = [];
      this.params.ProbabilisticConn.exclCorrBelow = 0; %i.e. set to 0.1 to exclude small corr values
      this.params.ProbabilisticConn.exclCorrWithPValueAbove = []; %i.e. set to 0.05
      this.params.ProbabilisticConn.FDRalpha = []; %i.e. set to 0.05 false-discovery-ratio
      this.params.ProbabilisticConn.maxdx = Inf;
      this.params.ProbabilisticConn.maxdy = Inf;
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      
      if ~isempty(this.params.ProbabilisticConn.exclCorrWithPValueAbove) || ~isempty(this.params.ProbabilisticConn.FDRalpha)
        connProb = load(fullfile(this.workpath,this.params.ProbabilisticConn.inCorrFile));
        pValues = connProb.p;
        connProb = connProb.patchCorrCov;
      else
        connProb = load(fullfile(this.workpath,this.params.ProbabilisticConn.inCorrFile),'patchCorrCov');
        connProb = connProb.patchCorrCov;
      end
      
      if this.params.ProbabilisticConn.maxdx*2+1 < size(connProb,1)
        cut = (size(connProb,1)-(this.params.ProbabilisticConn.maxdx*2+1))/2;
        connProb = connProb(cut+1:end-cut,:,:,:);
        if ~isempty(this.params.ProbabilisticConn.exclCorrWithPValueAbove) || ~isempty(this.params.ProbabilisticConn.FDRalpha)
          pValues = pValues(cut+1:end-cut,:,:,:);
        end
      end
      if this.params.ProbabilisticConn.maxdy*2+1 < size(connProb,2)
        cut = (size(connProb,2)-(this.params.ProbabilisticConn.maxdy*2+1))/2;
        connProb = connProb(:,cut+1:end-cut,:,:);
        if ~isempty(this.params.ProbabilisticConn.exclCorrWithPValueAbove) || ~isempty(this.params.ProbabilisticConn.FDRalpha)
          pValues = pValues(:,cut+1:end-cut,:,:);
        end
      end
      
      N = size(connProb,4);
      
      dx = cell(N,1);
      dy = cell(N,1);
      f0 = cell(N,1);
      f1 = cell(N,1);
      w = cell(N,1);
      dxInh = cell(N,1);
      dyInh = cell(N,1);
      f0Inh = cell(N,1);
      f1Inh = cell(N,1);
      wInh = cell(N,1);
      
      
      
      mindx = -(size(connProb,1)-1)/2;
      mindy = -(size(connProb,2)-1)/2;
      
      if ~isempty(this.params.ProbabilisticConn.exclCorrWithPValueAbove)
        connProb(pValues > this.params.ProbabilisticConn.exclCorrWithPValueAbove) = 0;
      end
        
      for f0id=1:N
        disp(num2str(f0id))
        
        if this.params.ProbabilisticConn.exclSelfConn
          connProb(-mindx+1,-mindy+1,f0id,f0id) = 0;
        end
        
        if ~isempty(this.params.ProbabilisticConn.FDRalpha)
          connProb(:,:,:,f0id) = fdr(connProb(:,:,:,f0id),pValues(:,:,:,f0id),this.params.ProbabilisticConn.FDRalpha);
        end
        
        if ~isempty(this.params.ProbabilisticConn.probRadiusFcn)
          probRad = feval(this.params.ProbabilisticConn.probRadiusFcn,sqrt(bsxfun(@plus,((1:size(connProb,1)).^2)',((1:size(connProb,2)).^2))));
          probRad = probRad/sum(probRad(:));
          connProb(:,:,:,f0id) = bsxfun(@times,connProb(:,:,:,f0id),probRad);
        end
        
        %% excitatory connections:
        desiredConn = squeeze(connProb(:,:,:,f0id));
        desiredConn(desiredConn<this.params.ProbabilisticConn.exclCorrBelow) = 0;
        desiredConn(isnan(desiredConn)) = 1e-5;
        
        [f1{f0id},dx{f0id},dy{f0id}] = this.addConns(desiredConn,this.params.ProbabilisticConn.numExc);
        dx{f0id} = dx{f0id}+mindx-1;
        dy{f0id} = dy{f0id}+mindy-1;
        
        f0{f0id}=f0id*ones(size(f1{f0id}));
        w{f0id}=ones(size(f1{f0id}));
        
        %% inhibitory connection:
        desiredConn = - squeeze(connProb(:,:,:,f0id));
        desiredConn(desiredConn<this.params.ProbabilisticConn.exclCorrBelow) = 0;
        desiredConn(isnan(desiredConn)) = 1e-5;
        
        [f1Inh{f0id},dxInh{f0id},dyInh{f0id}] = this.addConns(desiredConn,this.params.ProbabilisticConn.numInh);
        dxInh{f0id} = dxInh{f0id}+mindx-1;
        dyInh{f0id} = dyInh{f0id}+mindy-1;
        
        f0Inh{f0id} = f0id*ones(size(f1Inh{f0id}));
        wInh{f0id} = -ones(size(f1Inh{f0id}));
      end
      
      W.dx = cat(1,cell2mat(dx),cell2mat(dxInh));
      W.dy = cat(1,cell2mat(dy),cell2mat(dyInh));
      W.f0 = cat(1,cell2mat(f0),cell2mat(f0Inh));
      W.f1 = cat(1,cell2mat(f1),cell2mat(f1Inh));
      W.w = cat(1,cell2mat(w),cell2mat(wInh));
      
      savedir = fullfile(this.workpath,this.params.ProbabilisticConn.outWeightsFolder);
      if this.numJobs > 1
        savedir = fullfile(savedir,num2str(this.currJobid));
      end
      mkdir(savedir);
      save(fullfile(savedir,'weights.mat'),'W')
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function [f1,dx,dy] = addConns(this,connProb,numConns)
      
        if this.params.ProbabilisticConn.useDiscretesample
          if sum(connProb(:))==0
            connProb(:) = 1;
          end
          IX = discretesample(connProb(:)/sum(connProb(:)),numConns)';
          if this.params.ProbabilisticConn.exclDoubleConns
            IX = unique(IX);
            while length(IX)~=numConns
              IX = cat(1,IX,discretesample(connProb(:)/sum(connProb(:)),numConns-length(IX))');
              IX = unique(IX);
            end
          end
        else
          connProb = connProb - rand(size(connProb));
          [~,IX] = sort(-connProb(:));
        end
        [dx,dy,f1] = ind2sub(size(connProb),IX(1:numConns));
    end
    
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

