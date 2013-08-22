classdef ApplyWeights < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ApplyWeights(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ApplyWeights.inActFolder = 'inAct'; %relative to the workpath
      this.params.ApplyWeights.inActFilenames = 'act.*.mat';
      this.params.ApplyWeights.inPhaseFolder = 'inPhase'; %relative to the workpath
      this.params.ApplyWeights.inPhaseFilenames = 'phase.*.mat';
      this.params.ApplyWeights.fileid = [];
      this.params.ApplyWeights.outPhaseFolder = 'outPhase'; %relative to the workpath
      this.params.ApplyWeights.weightsFile = 'AEWeights/weights.mat'; %relative to the workpath

      % If the following parameters are set, they will overwrite the 
      % specifications in the weightsFile
      this.params.ApplyWeights.inputSubtract = [];
      this.params.ApplyWeights.inputScaling = [];
      this.params.ApplyWeights.inputSubsampling = [];
      this.params.ApplyWeights.shiftOutputdims = [];
      this.params.ApplyWeights.convType = []; %{valid},same,full
      this.params.ApplyWeights.W = []; %linear weights
      this.params.ApplyWeights.b = []; %bias
      this.params.ApplyWeights.actFcn = []; %function handle
      this.params.ApplyWeights.actFcn2 = []; %function handle
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});

    end
    
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      %list available activity files which match the given string
      inputfolder = fullfile(this.workpath,this.params.ActPhaseCombine.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.ActPhaseCombine.inActFilenames );

      inputfolderPhase = fullfile(this.workpath,this.params.ActPhaseCombine.inPhaseFolder);
      [ pathlistPhase, filelistPhase ] = dirrec( inputfolderPhase,this.params.ActPhaseCombine.inPhaseFilenames );

      savedir = fullfile(this.workpath,this.params.ApplyWeights.outActFolder);
      
      if isempty(this.params.ApplyWeights.weightsFile)
        conn = struct();
      else
        if isstruct(this.params.ApplyWeights.weightsFile)
          conn = this.params.ApplyWeights.weightsFile;
        elseif ischar(this.params.ApplyWeights.weightsFile)
          conn = load(fullfile(this.workpath,this.params.ApplyWeights.weightsFile));
        end
      end
      if ~isempty(this.params.ApplyWeights.inputSubtract)
        conn.inputSubtract = this.params.ApplyWeights.inputSubtract;
      end
      if ~isempty(this.params.ApplyWeights.inputScaling)
        conn.inputScaling = this.params.ApplyWeights.inputScaling;
      end
      if ~isempty(this.params.ApplyWeights.inputSubsampling)
        conn.inputSubsampling = this.params.ApplyWeights.inputSubsampling;
      end
      if ~isempty(this.params.ApplyWeights.shiftOutputdims)
        conn.shiftOutputdims = this.params.ApplyWeights.shiftOutputdims;
      end
      if ~isempty(this.params.ApplyWeights.convType)
        conn.convType = this.params.ApplyWeights.convType;
      end
      if ~isempty(this.params.ApplyWeights.W)
        conn.W = this.params.ApplyWeights.W;
      end
      if ~isempty(this.params.ApplyWeights.b)
        conn.b = this.params.ApplyWeights.b;
      end
      if ~isempty(this.params.ApplyWeights.actFcn)
        conn.actFcn = this.params.ApplyWeights.actFcn;
      end
      if ~isempty(this.params.ApplyWeights.actFcn2)
        conn.actFcn2 = this.params.ApplyWeights.actFcn2;
      end
      
      if ~isempty(this.params.ApplyWeights.fileid)
        filelist = filelist(this.params.ApplyWeights.fileid);
        pathlist = pathlist(this.params.ApplyWeights.fileid);
        filelistPhase = filelistPhase(this.params.ActPhaseCombine.fileid);
        pathlistPhase = pathlistPhase(this.params.ActPhaseCombine.fileid);
      end
      
      for fileid=1:length(filelist)
        disp(fullfile(pathlist{fileid},filelist{fileid}));
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        phase = load(fullfile(pathlistPhase{fileid},filelistPhase{fileid}));
        phase = phase.phase;
        if iscell(act)
          for cellid=1:numel(act)
            if ~isempty(this.params.ActPhaseCombine.inPhaseFolder)
              act{cellid} = ApplyWeights.applyWeights(act{cellid}, phase{cellid}, conn);
            end
          end
        else
          if ~isempty(this.params.ActPhaseCombine.inPhaseFolder)
            act = ApplyWeights.applyWeights(act, phase, conn);
          end
        end
        savepath = fullfile(savedir, pathlist{fileid}(length(inputfolder)+2:end));
        mkdir(savepath);
        save(fullfile(savepath,filelist{fileid}),'-v7.3','act');
      end
      
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
  methods(Static)
    
    function act = applyWeights( act, phase, conn)
      
      if isfield(conn,'inputSubtract') && ~isempty(conn.inputSubtract)
        act = bsxfun(@minus,act,conn.inputSubtract);
      end
      
      if isfield(conn,'inputScaling') && ~isempty(conn.inputScaling)
        act = bsxfun(@times,act,conn.inputScaling);
      end

      
      if isfield(conn,'inputSubsampling') && ~isempty(conn.inputSubsampling) && conn.inputSubsampling > 1
        
        act = act(1:end-mod(size(act,1),conn.inputSubsampling),1:end-mod(size(act,2),conn.inputSubsampling),:);
        
        act = reshape(act, [...
          conn.inputSubsampling ...
          size(act,1)/conn.inputSubsampling ...
          conn.inputSubsampling ...
          size(act,2)/conn.inputSubsampling ...
          size(act,3)]);
        act = permute(act,[2 4 1 3 5]);
        act = reshape(act,[...
          size(act,1) ...
          size(act,2) ...
          size(act,3)*size(act,4)*size(act,5)]);
        
        conn.W = reshape(conn.W,[...
          conn.inputSubsampling ...
          size(conn.W,1)/conn.inputSubsampling ...
          conn.inputSubsampling ...
          size(conn.W,2)/conn.inputSubsampling ...
          size(conn.W,3)...
          size(conn.W,4)...
          size(conn.W,5)...
          size(conn.W,6)]);
        conn.W = permute(conn.W,[2 4 1 3 5 6 7 8]);
        conn.W = reshape(conn.W,[...
          size(conn.W,1) ...
          size(conn.W,2) ...
          size(conn.W,3)*size(conn.W,4)*size(conn.W,5) ...
          size(conn.W,6)...
          size(conn.W,7)...
          size(conn.W,8)]);
      end
      
      if isfield(conn,'W') && ~isempty(conn.W)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Use only absolute values of conn.W:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        conn.W = abs(conn.W);
        act = abs(act) .* exp(1i * phase);
        
        
        if isfield(conn,'shiftOutputdims') && ~isempty(conn.shiftOutputdims) && conn.shiftOutputdims
          
          % TODO implement tiled conv here:
          
        else
          
          dimW = size(conn.W);
          
          if isfield(conn,'convType') &&  ~isempty(conn.convType) && strcmp(conn.convType,'full')
            
            % TODO
            
          elseif isfield(conn,'convType') &&  ~isempty(conn.convType) && strcmp(conn.convType,'same')
            
            actFull = zeros([size(act,1)+dimW(1)-1 size(act,2)+dimW(2)-1 size(act,3)]);
            dim1start = ceil(dimW(1)/2);
            dim2start = ceil(dimW(2)/2);
            actFull(dim1start:dim1start+size(act,1)-1, dim2start:dim2start+size(act,2)-1, :) = act;
            act = actFull;
            
            % change W: move dimensions 4 and 5 into the feature dimension
            conn.W = reshape(conn.W,[dimW(1) dimW(2) dimW(3) dimW(4)*dimW(5)*dimW(6)]);
            act = conv3d(conn.W,act);
            % change act: extract dx and dy dimensions from dimension 3 and put them in dimensions 1 and 2
            act = reshape(act,[size(act,1) size(act,2) dimW(4) dimW(5) size(act,3)/(dimW(4)*dimW(5))]);
            act = permute(act,[3 1 4 2 5]);
            act = reshape(act,[size(act,1)*size(act,2) size(act,3)*size(act,4) size(act,5)]);
            
          else %valid convtype

            if dimW(1)==size(act,1) && dimW(2)==size(act,2) %check if convolution not necessary

              conn.W = reshape(conn.W,[dimW(1)*dimW(2)*dimW(3) dimW(4)*dimW(5)*dimW(6)]);
              act = act(:)' * conn.W;
              act = reshape(act,dimW(4:6));

            else

              % change W: move dimensions 4 and 5 into the feature dimension
              conn.W = reshape(conn.W,[dimW(1) dimW(2) dimW(3) dimW(4)*dimW(5)*dimW(6)]);
              act = conv3d(conn.W,act);
              % change act: extract dx and dy dimensions from dimension 3 and put them in dimensions 1 and 2
              act = reshape(act,[size(act,1) size(act,2) dimW(4) dimW(5) size(act,3)/(dimW(4)*dimW(5))]);
              act = permute(act,[3 1 4 2 5]);
              act = reshape(act,[size(act,1)*size(act,2) size(act,3)*size(act,4) size(act,5)]);

            end
            
          end
        end
      end
      
      if isfield(conn,'b') && ~isempty(conn.b)
        act = bsxfun(@plus,act,conn.b);
      end
      
      if isfield(conn,'actFcn') && ~isempty(conn.actFcn)
        act = feval(conn.actFcn,act);
      end

      if isfield(conn,'actFcn2') && ~isempty(conn.actFcn2)
        act = feval(conn.actFcn2,act);
      end
      
    end
    
  end
  
end

