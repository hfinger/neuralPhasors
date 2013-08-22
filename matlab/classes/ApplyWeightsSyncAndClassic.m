classdef ApplyWeightsSyncAndClassic < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    stats = struct()
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ApplyWeightsSyncAndClassic(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ApplyWeightsSyncAndClassic.inActFolder = 'inAct'; %relative to the workpath
      this.params.ApplyWeightsSyncAndClassic.inActFilenames = 'act.*.mat';
      this.params.ApplyWeightsSyncAndClassic.inPhaseFolder = 'inPhase'; %relative to the workpath
      this.params.ApplyWeightsSyncAndClassic.inPhaseFilenames = 'phase.*.mat';
      this.params.ApplyWeightsSyncAndClassic.fileid = [];
      this.params.ApplyWeightsSyncAndClassic.outActFolder = 'outAct'; %relative to the workpath
      this.params.ApplyWeightsSyncAndClassic.outPhaseFolder = 'outPhase'; %relative to the workpath
      this.params.ApplyWeightsSyncAndClassic.outStatsFolder = []; %relative to the workpath
      this.params.ApplyWeightsSyncAndClassic.weightsFile = 'AEWeights/weights.mat'; %relative to the workpath
      this.params.ApplyWeightsSyncAndClassic.plotPdf = false; %if true then plot a pdf
      this.params.ApplyWeightsSyncAndClassic.weightSyncTerm = 0.5; % between 0 (classic) and 1 (sync) to weight between classic and sync term
      this.params.ApplyWeightsSyncAndClassic.useAbsWeight = false;
      this.params.ApplyWeightsSyncAndClassic.calcStatsBins = 0;

      % If the following parameters are set, they will overwrite the 
      % specifications in the weightsFile
      this.params.ApplyWeightsSyncAndClassic.inputSubsampling = [];
      this.params.ApplyWeightsSyncAndClassic.shiftOutputdims = [];
      this.params.ApplyWeightsSyncAndClassic.convType = []; %{valid},same,full
      this.params.ApplyWeightsSyncAndClassic.W = []; %linear weights
      this.params.ApplyWeightsSyncAndClassic.b = []; % only applied to output act and not to output phase
      this.params.ApplyWeightsSyncAndClassic.actFcn = []; % only applied to output act and not to output phase
      this.params.ApplyWeightsSyncAndClassic.actFcn2 = []; % only applied to output act and not to output phase
      this.params.ApplyWeightsSyncAndClassic.outputUpSampling = [];
      
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
      inputfolder = fullfile(this.workpath,this.params.ApplyWeightsSyncAndClassic.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.ApplyWeightsSyncAndClassic.inActFilenames );
      
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.inPhaseFolder)
        inputfolderPhase = fullfile(this.workpath,this.params.ApplyWeightsSyncAndClassic.inPhaseFolder);
        [ pathlistPhase, filelistPhase ] = dirrec( inputfolderPhase,this.params.ApplyWeightsSyncAndClassic.inPhaseFilenames );
      end
      
      if isempty(this.params.ApplyWeightsSyncAndClassic.weightsFile)
        conn = struct();
      else
        if isstruct(this.params.ApplyWeightsSyncAndClassic.weightsFile)
          conn = this.params.ApplyWeightsSyncAndClassic.weightsFile;
        elseif ischar(this.params.ApplyWeightsSyncAndClassic.weightsFile)
          conn = load(fullfile(this.workpath,this.params.ApplyWeightsSyncAndClassic.weightsFile));
        end
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.inputSubsampling)
        conn.inputSubsampling = this.params.ApplyWeightsSyncAndClassic.inputSubsampling;
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.shiftOutputdims)
        conn.shiftOutputdims = this.params.ApplyWeightsSyncAndClassic.shiftOutputdims;
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.convType)
        conn.convType = this.params.ApplyWeightsSyncAndClassic.convType;
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.W)
        if this.params.ApplyWeightsSyncAndClassic.W==false
          conn.W = [];
        else
          conn.W = this.params.ApplyWeightsSyncAndClassic.W;
        end
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.b)
        conn.b = this.params.ApplyWeightsSyncAndClassic.b;
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.actFcn)
        conn.actFcn = this.params.ApplyWeightsSyncAndClassic.actFcn;
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.actFcn2)
        conn.actFcn2 = this.params.ApplyWeightsSyncAndClassic.actFcn2;
      end
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.outputUpSampling)
        conn.outputUpSampling = this.params.ApplyWeightsSyncAndClassic.outputUpSampling;
      end
      
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.fileid)
        filelist = filelist(this.params.ApplyWeightsSyncAndClassic.fileid);
        pathlist = pathlist(this.params.ApplyWeightsSyncAndClassic.fileid);
        if ~isempty(this.params.ApplyWeightsSyncAndClassic.inPhaseFolder)
          filelistPhase = filelistPhase(this.params.ApplyWeightsSyncAndClassic.fileid);
          pathlistPhase = pathlistPhase(this.params.ApplyWeightsSyncAndClassic.fileid);
        end
      end
      
      if this.params.ApplyWeightsSyncAndClassic.useAbsWeight
        conn.W = abs(conn.W);
      end
      
      for fileid=1:length(filelist)
        disp(['fileid: ' num2str(fileid) ' filename: ' fullfile(pathlist{fileid},filelist{fileid})]);
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        if ~isempty(this.params.ApplyWeightsSyncAndClassic.inPhaseFolder)
          phase = load(fullfile(pathlistPhase{fileid},filelistPhase{fileid}));
          phase = phase.phase;
        else
          phase = [];
        end
        if iscell(act)
          for cellid=1:numel(act)
            [this, act{cellid}, phase{cellid}] = this.applyWeights(act{cellid}, phase, conn);
          end
        else
          [this, act, phase] = this.applyWeights(act, phase, conn);
        end
        
        
        
        if ~isempty(this.params.ApplyWeightsSyncAndClassic.outPhaseFolder)
          savedirPhase = fullfile(this.workpath,this.params.ApplyWeightsSyncAndClassic.outPhaseFolder);
          if this.numJobs > 1
            savedirPhase = fullfile(savedirPhase,num2str(this.currJobid));
          end
          savepathPhase = fullfile(savedirPhase, pathlistPhase{fileid}(length(inputfolderPhase)+2:end));
          mkdir(savepathPhase);
          save(fullfile(savepathPhase,filelistPhase{fileid}),'-v7.3','phase');
        end
        
        if ~isempty(this.params.ApplyWeightsSyncAndClassic.outActFolder)
          savepath = fullfile(this.workpath,this.params.ApplyWeightsSyncAndClassic.outActFolder);
          if this.numJobs > 1
            savepath = fullfile(savepath,num2str(this.currJobid));
          end
          savepath = fullfile(savepath, pathlist{fileid}(length(inputfolder)+2:end));
      
          mkdir(savepath);
          save(fullfile(savepath,filelist{fileid}),'-v7.3','act');
          
          if this.params.ApplyWeightsSyncAndClassic.plotPdf
            fDim1 = floor(sqrt(size(act,3)));
            fDim2 = round(size(act,3)/fDim1);
            while size(act,3) ~= fDim1 * fDim2
              fDim1 = fDim1 + 1;
              fDim2 = round(size(act,3)/fDim1);
            end
            pdfName = fullfile(savepath,filelist{fileid});
            if strcmp(pdfName(end-3:end),'.mat')
              pdfName(end-3:end) = [];
            end
            plotFeatures( reshape(act,[size(act,1) size(act,2) fDim1 fDim2]), pdfName, 'jet', 5, [], [], false, true, false, false );
          end
        end
        
      end
      
      if ~isempty(this.params.ApplyWeightsSyncAndClassic.outStatsFolder)
        savepathStats = fullfile(this.workpath,this.params.ApplyWeightsSyncAndClassic.outStatsFolder);
        if this.numJobs > 1
          savepathStats = fullfile(savepathStats,num2str(this.currJobid));
        end
        mkdir(savepathStats);
        stats = this.stats;
        save(fullfile(savepathStats,'stats.mat'),'-v7.3','stats');
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
    
    
    function [this, act, phaseOut] = applyWeights( this, act, phase, conn)
      
      
      
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
      
      if isfield(conn,'outputUpSampling') && ~isempty(conn.outputUpSampling) && conn.outputUpSampling > 1
        conn.W = reshape(conn.W,[...
          size(conn.W,1) ...
          size(conn.W,2) ...
          size(conn.W,3) ...
          conn.outputUpSampling...
          size(conn.W,4)/conn.outputUpSampling...
          conn.outputUpSampling...
          size(conn.W,5)/conn.outputUpSampling...
          size(conn.W,6)]);
        conn.W = permute(conn.W,[1 2 3 5 7 4 6 8]);
        conn.W = reshape(conn.W,[...
          size(conn.W,1) ...
          size(conn.W,2) ...
          size(conn.W,3) ...
          size(conn.W,4)...
          size(conn.W,5)...
          size(conn.W,6)*size(conn.W,7)*size(conn.W,8)]);
      end
      
      %% sync term:
      actSync = ApplyWeightsSyncAndClassic.applyComplexConv(conn, act .* exp(1i * phase) );
      phaseOut = angle(actSync);
      idsZeros = find(actSync(:)==0);
      phaseOut(idsZeros) = 2*pi*rand(length(idsZeros),1);
      
      %% calc statistics
      if this.params.ApplyWeightsSyncAndClassic.calcStatsBins > 0
        assert(strcmp(conn.convType,'valid'))
        
        if ~isfield(this.stats,'allConns')
          this.stats.allConns.accumRelConnStrength = zeros(this.params.ApplyWeightsSyncAndClassic.calcStatsBins,size(actSync,3));
          this.stats.allConns.accumConnStrength = zeros(this.params.ApplyWeightsSyncAndClassic.calcStatsBins,size(actSync,3));
          this.stats.allConns.accumAct = zeros(this.params.ApplyWeightsSyncAndClassic.calcStatsBins,size(actSync,3));
          this.stats.allConns.accumNumNeurons = zeros(this.params.ApplyWeightsSyncAndClassic.calcStatsBins,size(actSync,3));
          this.stats.allConns.accumWeights = zeros(this.params.ApplyWeightsSyncAndClassic.calcStatsBins,size(actSync,3));
          this.stats.posConn = this.stats.allConns;
          this.stats.negConn = this.stats.allConns;
        end
        
        accumDims = [this.params.ApplyWeightsSyncAndClassic.calcStatsBins 1];
        
        for f=1:size(actSync,3)
          disp(['f=' num2str(f)]);
          
          currConnW = conn.W(:,:,:,1,1,f);
          posConns = (currConnW(:)>0);
          negConns = (currConnW(:)<0);
          
          numRand = 100;
          xRand = randi(size(actSync,1),1,numRand);
          yRand = randi(size(actSync,2),1,numRand);
          for j=1:numRand
            x=xRand(j);
            y=yRand(j);
          %for x=1:size(actSync,1)
            %disp(['x=' num2str(x)]);
            %for y=1:size(actSync,2)
            
              currAct = act(x:x+size(conn.W,1)-1,y:y+size(conn.W,2)-1,:);
              currPhase = phase(x:x+size(conn.W,1)-1,y:y+size(conn.W,2)-1,:);
              currPhaseDiff = mod(bsxfun(@minus,currPhase,phaseOut(x,y,f)),2*pi);
              currPhaseBinId = 1+floor(currPhaseDiff*this.params.ApplyWeightsSyncAndClassic.calcStatsBins / (2*pi));
              currConnStrength = currConnW(:) .* currAct(:);
              
              this.stats.allConns.accumRelConnStrength(:,f) = this.stats.allConns.accumRelConnStrength(:,f) + accumarray(currPhaseBinId(:), currConnStrength(:).*cos(currPhaseDiff(:)) , accumDims);
              this.stats.allConns.accumConnStrength(:,f) = this.stats.allConns.accumConnStrength(:,f) + accumarray(currPhaseBinId(:), currConnStrength(:), accumDims);
              this.stats.allConns.accumAct(:,f) = this.stats.allConns.accumAct(:,f) + accumarray(currPhaseBinId(:), currAct(:), accumDims);
              this.stats.allConns.accumNumNeurons(:,f) = this.stats.allConns.accumNumNeurons(:,f) + accumarray(currPhaseBinId(:), ones(size(currConnStrength(:))), accumDims);              
              this.stats.allConns.accumWeights(:,f) = this.stats.allConns.accumWeights(:,f) + accumarray(currPhaseBinId(:), currConnW(:), accumDims);              
              
              this.stats.posConn.accumRelConnStrength(:,f) = this.stats.posConn.accumRelConnStrength(:,f) + accumarray(currPhaseBinId(posConns), currConnStrength(posConns).*cos(currPhaseDiff(posConns)) , accumDims);
              this.stats.posConn.accumConnStrength(:,f) = this.stats.posConn.accumConnStrength(:,f) + accumarray(currPhaseBinId(posConns), currConnStrength(posConns), accumDims);
              this.stats.posConn.accumAct(:,f) = this.stats.posConn.accumAct(:,f) + accumarray(currPhaseBinId(posConns), currAct(posConns), accumDims);
              this.stats.posConn.accumNumNeurons(:,f) = this.stats.posConn.accumNumNeurons(:,f) + accumarray(currPhaseBinId(posConns), ones(size(currConnStrength(posConns))), accumDims);              
              this.stats.posConn.accumWeights(:,f) = this.stats.posConn.accumWeights(:,f) + accumarray(currPhaseBinId(posConns), currConnW(posConns), accumDims);              
              
              this.stats.negConn.accumRelConnStrength(:,f) = this.stats.negConn.accumRelConnStrength(:,f) + accumarray(currPhaseBinId(negConns), currConnStrength(negConns).*cos(currPhaseDiff(negConns)) , accumDims);
              this.stats.negConn.accumConnStrength(:,f) = this.stats.negConn.accumConnStrength(:,f) + accumarray(currPhaseBinId(negConns), currConnStrength(negConns), accumDims);
              this.stats.negConn.accumAct(:,f) = this.stats.negConn.accumAct(:,f) + accumarray(currPhaseBinId(negConns), currAct(negConns), accumDims);
              this.stats.negConn.accumNumNeurons(:,f) = this.stats.negConn.accumNumNeurons(:,f) + accumarray(currPhaseBinId(negConns), ones(size(currConnStrength(negConns))), accumDims);              
              this.stats.negConn.accumWeights(:,f) = this.stats.negConn.accumWeights(:,f) + accumarray(currPhaseBinId(negConns), currConnW(negConns), accumDims);              
              
            %end
          %end
          end
        end
        
        
        
      end
      
      %% combine sync term with classic term:
      alpha = this.params.ApplyWeightsSyncAndClassic.weightSyncTerm;
      act = alpha * abs(actSync) + (1-alpha) * ApplyWeightsSyncAndClassic.applyComplexConv(conn, act);
      
      
      
      assert(isreal(act)); % this is only the rate variable and not phase!
          
      if isfield(conn,'b') && ~isempty(conn.b)
        act = bsxfun(@plus,act,conn.b);
      end
      
      if isfield(conn,'actFcn') && ~isempty(conn.actFcn)
        act = feval(conn.actFcn,act);
      end

      if isfield(conn,'actFcn2') && ~isempty(conn.actFcn2)
        act = feval(conn.actFcn2,act);
      end
      
      if isfield(conn,'outputUpSampling') && ~isempty(conn.outputUpSampling) && conn.outputUpSampling > 1
        act = reshape(act, [...
          size(act,1) ...
          size(act,2) ...
          conn.outputUpSampling ...
          conn.outputUpSampling ...
          size(act,3)/(conn.outputUpSampling^2)]);
        act = permute(act,[3 1 4 2 5]);
        act = reshape(act,[...
          size(act,1)*size(act,2) ...
          size(act,3)*size(act,4) ...
          size(act,5)]);
        
        phaseOut = reshape(phaseOut, [...
          size(phaseOut,1) ...
          size(phaseOut,2) ...
          conn.outputUpSampling ...
          conn.outputUpSampling ...
          size(phaseOut,3)/(conn.outputUpSampling^2)]);
        phaseOut = permute(phaseOut,[3 1 4 2 5]);
        phaseOut = reshape(phaseOut,[...
          size(phaseOut,1)*size(phaseOut,2) ...
          size(phaseOut,3)*size(phaseOut,4) ...
          size(phaseOut,5)]);
        
        
      end
      
    end
    
    
  end
  methods(Static)
    
    
    
    %% act can be complex
    function act = applyComplexConv(conn, act)
      if isfield(conn,'W') && ~isempty(conn.W)
        if isfield(conn,'shiftOutputdims') && ~isempty(conn.shiftOutputdims) && conn.shiftOutputdims
          
          % TODO implement tiled conv here:
          
          tileSize1=size(conn.W,1)-1;
          tileSize2=size(conn.W,2)-1;
          fIn = size(conn.W,3);
          fOut = size(conn.W,6);
          this.tiledConvF = TiledConv(tileSize1,tileSize2,fIn,fOut);
        
          act = reshape(act,[size(act,1) size(act,2) tileSize1 tileSize2 fIn]);
          act = permute(act,[3 1 4 2 5]);
          act = reshape(act,[size(act,1)*size(act,2) size(act,3)*size(act,4) size(act,5)]);
          
          cutSize1 = mod(size(act,1),tileSize1);
          cutSize2 = mod(size(act,2),tileSize2);
          
          act = this.tiledConvF.convW(conn.W,act(1:end-cutSize1,1:end-cutSize2,:),'valid');
          act = permute(act,[3 1 4 2 5]);
          act = reshape(act,[size(act,1)*size(act,2) size(act,3)*size(act,4) size(act,5)]);
          
          
          
          
        else
          
          dimW = size(conn.W);
          dimW(end+1:6) = 1;
          
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
      
    end
    
  end
  
end

