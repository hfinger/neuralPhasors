classdef CombineWeights < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = CombineWeights(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.CombineWeights.weightsFileIn1WithoutActFcn = 'WhiteningWeights/weights.mat'; %relative to the workpath
      this.params.CombineWeights.weightsFileIn2WithActFCN = 'AEWeights/weights.mat'; %relative to the workpath
      this.params.CombineWeights.weightsFileOutCombined = 'weightsCombined'; %relative to the workpath
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      
      conn1 = load(fullfile(this.workpath,this.params.CombineWeights.weightsFileIn1WithoutActFcn));
      conn2 = load(fullfile(this.workpath,this.params.CombineWeights.weightsFileIn2WithActFCN));
      connCombined = struct();
      
      %% from conn1:
      if isfield(conn1,'inputSubtract') && ~isempty(conn1.inputSubtract)
        connCombined.inputSubtract = conn1.inputSubtract;
      end
      if isfield(conn1,'inputScaling') && ~isempty(conn1.inputScaling)
        connCombined.inputScaling = conn1.inputScaling;
      end
      if isfield(conn1,'inputSubsampling') && ~isempty(conn1.inputSubsampling)
        connCombined.inputSubsampling = conn1.inputSubsampling;
      end
      
      %% combine conn1 and conn2:
      assert(size(conn1.W,4)==1 && size(conn1.W,5)==1)
      assert(size(conn2.W,4)==1 && size(conn2.W,5)==1)
      assert(size(conn1.W,6)==size(conn2.W,3)) % center dim must be the same
      newWdim1 = 1 + size(conn2.W,1)-1 + size(conn1.W,1)-1;
      newWdim2 = 1 + size(conn2.W,2)-1 + size(conn1.W,2)-1;
      newWdim3 = size(conn1.W,3);
      newWdim4 = 1;
      newWdim5 = 1;
      newWdim6 = size(conn2.W,6);
      connCombined.W = zeros(newWdim1,newWdim2,newWdim3,newWdim4,newWdim5,newWdim6);
      for indim=1:newWdim3
        for outdim=1:newWdim6
          for centerdim=1:size(conn1.W,6)
            connCombined.W(:,:,indim,1,1,outdim) = connCombined.W(:,:,indim,1,1,outdim) + ...
              conv2( conn1.W(end:-1:1,end:-1:1,indim,1,1,centerdim), conn2.W(:,:,centerdim,1,1,outdim) );
          end
        end
      end
      
      %% from conn2:
      if isfield(conn2,'shiftOutputdims') && ~isempty(conn2.shiftOutputdims)
        connCombined.shiftOutputdims = conn2.shiftOutputdims;
      end
      if isfield(conn2,'b') && ~isempty(conn2.b)
        connCombined.b = conn2.b;
      end
      if isfield(conn2,'actFcn') && ~isempty(conn2.actFcn)
        connCombined.actFcn = conn2.actFcn;
      end
      if isfield(conn2,'actFcn2') && ~isempty(conn2.actFcn2)
        connCombined.actFcn2 = conn2.actFcn2;
      end
      
      savepath = fullfile(this.workpath,this.params.CombineWeights.weightsFileOutCombined);
      if this.numJobs > 1
        savepath = fullfile(savepath,num2str(this.currJobid));
      end
      mkdir(savepath);
      
      save(fullfile(savepath,'forwConn.mat'),'-struct','connCombined')
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
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

