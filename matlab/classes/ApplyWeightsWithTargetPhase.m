classdef ApplyWeightsWithTargetPhase < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ApplyWeightsWithTargetPhase(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ApplyWeightsWithTargetPhase.inActFolder = 'inAct'; %relative to the workpath
      this.params.ApplyWeightsWithTargetPhase.inActFilenames = 'act.*.mat';
      this.params.ApplyWeightsWithTargetPhase.inPhaseFolder = []; %relative to the workpath
      this.params.ApplyWeightsWithTargetPhase.inPhaseFilenames = [];
      this.params.ApplyWeightsWithTargetPhase.inTargetPhaseFolder = []; %relative to the workpath
      this.params.ApplyWeightsWithTargetPhase.inTargetPhaseFilenames = [];
      this.params.ApplyWeightsWithTargetPhase.fileid = [];
      this.params.ApplyWeightsWithTargetPhase.outActFolder = 'outAct'; %relative to the workpath
      this.params.ApplyWeightsWithTargetPhase.weightsFile = 'AEWeights/weights.mat'; %relative to the workpath

      % If the following parameters are set, they will overwrite the 
      % specifications in the weightsFile
      this.params.ApplyWeightsWithTargetPhase.inputSubtract = [];
      this.params.ApplyWeightsWithTargetPhase.inputScaling = [];
      this.params.ApplyWeightsWithTargetPhase.inputSubsampling = [];
      this.params.ApplyWeightsWithTargetPhase.shiftOutputdims = [];
      this.params.ApplyWeightsWithTargetPhase.convType = []; %same
      this.params.ApplyWeightsWithTargetPhase.W = []; %linear weights
      this.params.ApplyWeightsWithTargetPhase.b = []; %bias
      this.params.ApplyWeightsWithTargetPhase.actFcn = []; %function handle
      this.params.ApplyWeightsWithTargetPhase.actFcn2 = []; %function handle
      
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
      inputfolder = fullfile(this.workpath,this.params.ApplyWeightsWithTargetPhase.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.ApplyWeightsWithTargetPhase.inActFilenames );
      inputfolderPhase = fullfile(this.workpath,this.params.ApplyWeightsWithTargetPhase.inPhaseFolder);
      [ pathlistPhase, filelistPhase ] = dirrec( inputfolderPhase,this.params.ApplyWeightsWithTargetPhase.inPhaseFilenames );
      inputfolderTargetPhase = fullfile(this.workpath,this.params.ApplyWeightsWithTargetPhase.inTargetPhaseFolder);
      [ pathlistTargetPhase, filelistTargetPhase ] = dirrec( inputfolderPhase,this.params.ApplyWeightsWithTargetPhase.inTargetPhaseFilenames );
      
      if isempty(this.params.ApplyWeightsWithTargetPhase.weightsFile)
        conn = struct();
      else
        if isstruct(this.params.ApplyWeightsWithTargetPhase.weightsFile)
          conn = this.params.ApplyWeightsWithTargetPhase.weightsFile;
        elseif ischar(this.params.ApplyWeightsWithTargetPhase.weightsFile)
          conn = load(fullfile(this.workpath,this.params.ApplyWeightsWithTargetPhase.weightsFile));
        end
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.inputSubtract)
        conn.inputSubtract = this.params.ApplyWeightsWithTargetPhase.inputSubtract;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.inputScaling)
        conn.inputScaling = this.params.ApplyWeightsWithTargetPhase.inputScaling;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.inputSubsampling)
        conn.inputSubsampling = this.params.ApplyWeightsWithTargetPhase.inputSubsampling;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.shiftOutputdims)
        conn.shiftOutputdims = this.params.ApplyWeightsWithTargetPhase.shiftOutputdims;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.convType)
        conn.convType = this.params.ApplyWeightsWithTargetPhase.convType;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.W)
        conn.W = this.params.ApplyWeightsWithTargetPhase.W;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.b)
        conn.b = this.params.ApplyWeightsWithTargetPhase.b;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.actFcn)
        conn.actFcn = this.params.ApplyWeightsWithTargetPhase.actFcn;
      end
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.actFcn2)
        conn.actFcn2 = this.params.ApplyWeightsWithTargetPhase.actFcn2;
      end
      
      if ~isempty(this.params.ApplyWeightsWithTargetPhase.fileid)
        filelist = filelist(this.params.ApplyWeightsWithTargetPhase.fileid);
        pathlist = pathlist(this.params.ApplyWeightsWithTargetPhase.fileid);
        filelistPhase = filelistPhase(this.params.ApplyWeightsWithTargetPhase.fileid);
        pathlistPhase = pathlistPhase(this.params.ApplyWeightsWithTargetPhase.fileid);
        filelistTargetPhase = filelistTargetPhase(this.params.ApplyWeightsWithTargetPhase.fileid);
        pathlistTargetPhase = pathlistTargetPhase(this.params.ApplyWeightsWithTargetPhase.fileid);
      end
      
      for fileid=1:length(filelist)
        disp(fullfile(pathlist{fileid},filelist{fileid}));
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        phase = load(fullfile(pathlistPhase{fileid},filelistPhase{fileid}));
        phase = phase.phase;
        targetPhase = load(fullfile(pathlistTargetPhase{fileid},filelistTargetPhase{fileid}));
        targetPhase = targetPhase.phase;
        if iscell(act)
          for cellid=1:numel(act)
            act{cellid} = ApplyWeightsWithTargetPhase.applyWeights(act{cellid}, phase, conn, targetPhase);
          end
        else
          act = ApplyWeightsWithTargetPhase.applyWeights(act, phase, conn, targetPhase);
        end
        
        if ~isempty(this.params.ApplyWeightsWithTargetPhase.outActFolder)
          savedir = fullfile(this.workpath,this.params.ApplyWeightsWithTargetPhase.outActFolder);
          savepath = fullfile(savedir, pathlist{fileid}(length(inputfolder)+2:end));
          mkdir(savepath);
          save(fullfile(savepath,filelist{fileid}),'-v7.3','act');
        end
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
    
    function act = applyWeights( act, phase, conn, targetPhase)
      
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
        
        phase = phase(1:end-mod(size(phase,1),conn.inputSubsampling),1:end-mod(size(phase,2),conn.inputSubsampling),:);
        phase = reshape(phase, [...
          conn.inputSubsampling ...
          size(phase,1)/conn.inputSubsampling ...
          conn.inputSubsampling ...
          size(phase,2)/conn.inputSubsampling ...
          size(phase,3)]);
        phase = permute(phase,[2 4 1 3 5]);
        phase = reshape(phase,[...
          size(phase,1) ...
          size(phase,2) ...
          size(phase,3)*size(phase,4)*size(phase,5)]);
        
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
        if isfield(conn,'shiftOutputdims') && ~isempty(conn.shiftOutputdims) && conn.shiftOutputdims
          
          % TODO implement tiled conv here:
          
        else
          
          dimW = size(conn.W);
          
          if isfield(conn,'convType') &&  ~isempty(conn.convType) && strcmp(conn.convType,'full')
            
            % TODO
            
          elseif isfield(conn,'convType') &&  ~isempty(conn.convType) && strcmp(conn.convType,'same')
            
            %implement convolution: z2 = W*data*cos(dataPhase-targetPhase)
            %use trigonometric identity: cos(x-y) = cos(x)cos(y) + sin(x)sin(y)
            %implement convolution: z2 = W*data*cos(dataPhase)*cos(targetPhase) + W*data*sin(dataPhase)*sin(targetPhase)

            dim1start = ceil(dimW(1)/2);
            dim2start = ceil(dimW(2)/2);
            
            % change W: move dimensions 4 and 5 into the feature dimension
            conn.W = reshape(conn.W,[dimW(1) dimW(2) dimW(3) dimW(4)*dimW(5)*dimW(6)]);
            
            actFullC = zeros([size(act,1)+dimW(1)-1 size(act,2)+dimW(2)-1 size(act,3)]);
            actFullS = zeros([size(act,1)+dimW(1)-1 size(act,2)+dimW(2)-1 size(act,3)]);
            
            actFullC(dim1start:dim1start+size(act,1)-1, dim2start:dim2start+size(act,2)-1, :) = act.*cos(phase);
            actFullS(dim1start:dim1start+size(act,1)-1, dim2start:dim2start+size(act,2)-1, :) = act.*sin(phase);

            for i=1:2
              if i==1
                act = conv3d(conn.W,actFullC);
              else
                act = conv3d(conn.W,actFullS);
              end
              % change act: extract dx and dy dimensions from dimension 3 and put them in dimensions 1 and 2
              act = reshape(act,[size(act,1) size(act,2) dimW(4) dimW(5) size(act,3)/(dimW(4)*dimW(5))]);
              act = permute(act,[3 1 4 2 5]);
              act = reshape(act,[size(act,1)*size(act,2) size(act,3)*size(act,4) size(act,5)]);
              if i==1
                actOut = act.*cos(targetPhase);
              else
                actOut = actOut + act.*sin(targetPhase);
              end
            end
            act = actOut;
            
          else %valid convtype

            %% TODO: implement with targetPhase
%             if dimW(1)==size(act,1) && dimW(2)==size(act,2) %check if convolution not necessary
% 
%               conn.W = reshape(conn.W,[dimW(1)*dimW(2)*dimW(3) dimW(4)*dimW(5)*dimW(6)]);
%               act = act(:)' * conn.W;
%               act = reshape(act,dimW(4:6));
% 
%             else
% 
%               % change W: move dimensions 4 and 5 into the feature dimension
%               conn.W = reshape(conn.W,[dimW(1) dimW(2) dimW(3) dimW(4)*dimW(5)*dimW(6)]);
%               act = conv3d(conn.W,act);
%               % change act: extract dx and dy dimensions from dimension 3 and put them in dimensions 1 and 2
%               act = reshape(act,[size(act,1) size(act,2) dimW(4) dimW(5) size(act,3)/(dimW(4)*dimW(5))]);
%               act = permute(act,[3 1 4 2 5]);
%               act = reshape(act,[size(act,1)*size(act,2) size(act,3)*size(act,4) size(act,5)]);
% 
%             end
            
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
    
    function z2 = conv3dTargetPhase(W,data,dataPhase,targetPhase)
      %data has dimensions x,y,fin
      %W has dimensions dx,dy,fin,fout
      %z2 has dimensions x,y,fout
      %dataPhase has dimensions x,y,fin
      %targetPhase has dimensions x,y,fout
      
      %implement as convolution: z2 = W*data*cos(dataPhase-targetPhase)
      %use trigonometric identity: cos(x-y) = cos(x)cos(y) + sin(x)sin(y)
      %implement as convolution: z2 = W*data*cos(dataPhase)*cos(targetPhase) + W*data*sin(dataPhase)*sin(targetPhase)
      
      z2 = zeros([size(data,1)-size(W,1)+1 size(data,2)-size(W,2)+1 size(W,4)]); %this is initially the cos part
      z2sin = zeros([size(data,1)-size(W,1)+1 size(data,2)-size(W,2)+1 size(W,4)]);
      for f1=1:size(data,3)
        blockFilters = reshape(W(:,:,f1,:),[size(W,1)*size(W,2) size(W,4)]);
        
        %add cos part:
        blockImg = im2col(data(:,:,f1).*cos(dataPhase(:,:,f1)),[size(W,1) size(W,2)]);
        z2cos = z2cos + reshape(blockImg'*blockFilters,[size(data,1)-size(W,1)+1 size(data,2)-size(W,2)+1 size(W,4)]);
        
        %add sin part:
        blockImg = im2col(data(:,:,f1).*sin(dataPhase(:,:,f1)),[size(W,1) size(W,2)]);
        z2sin = z2sin + reshape(blockImg'*blockFilters,[size(data,1)-size(W,1)+1 size(data,2)-size(W,2)+1 size(W,4)]);
      end
      
      z2 = z2.*cos(targetPhase);
      z2 = z2 + z2sin.*sin(targetPhase);
    end
    
  end
  
end

