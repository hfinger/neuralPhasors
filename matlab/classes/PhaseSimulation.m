classdef PhaseSimulation < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
    tempsavepath
    savepath
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = PhaseSimulation(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.PhaseSimulation.inActFolder = 'labelMeWhite'; %relative to the workpath
      this.params.PhaseSimulation.inActFilenames = 'act.*.mat';
      this.params.PhaseSimulation.inPhaseFolder = []; %relative to the workpath (should correspond to inActFolder)
      this.params.PhaseSimulation.inPhaseFilenames = []; %relative to the workpath (should correspond to inActFilenames)
      this.params.PhaseSimulation.inFileid = 1;
      this.params.PhaseSimulation.inCellid = 1;
      this.params.PhaseSimulation.inConnFilename = 'conn/weights.mat'; %relative to the workpath
      this.params.PhaseSimulation.inPhaseFilename = []; %relative to the workpath (if init all simulations with same phase)
      this.params.PhaseSimulation.outPhaseFolder = 'labelMePhase'; %relative to the workpath
      
      this.params.PhaseSimulation.noiseLevel = 0;
      this.params.PhaseSimulation.noiseEMAconst = 0;
      this.params.PhaseSimulation.tmax = 20;
      this.params.PhaseSimulation.dt = 1;
      this.params.PhaseSimulation.fixedPhaseDelay = 0;
      this.params.PhaseSimulation.odeSolver = 'ode1'; %ode1,ode2,ode3,ode4,ode23
      this.params.PhaseSimulation.weightAll = 1;
      this.params.PhaseSimulation.weightInh = 1;
      this.params.PhaseSimulation.weightExc = 1;
      this.params.PhaseSimulation.saveintervalPhase = 1;
      this.params.PhaseSimulation.saveintervalMeanPhase = Inf;
      this.params.PhaseSimulation.saveintervalMeanWeightedPhase = Inf;
      this.params.PhaseSimulation.saveintervalMeanWeightedPhasePng = Inf;
      this.params.PhaseSimulation.plotPhase = false;
      this.params.PhaseSimulation.maxdphase = 0.5;
      this.params.PhaseSimulation.weightInhRadialFcn = []; % i.e. @(r) min(0,5-r)
      this.params.PhaseSimulation.weightGlobInh = 0;
      this.params.PhaseSimulation.sampleActOnce = false;
      this.params.PhaseSimulation.sampleActRepeated = false;
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      if ischar(this.params.PhaseSimulation.inActFolder)
        inputfolder = fullfile(this.workpath,this.params.PhaseSimulation.inActFolder);
        [ pathlist, filelist ] = dirrec( inputfolder,this.params.PhaseSimulation.inActFilenames );
        if ~isempty(this.params.PhaseSimulation.inFileid)
          filelist = filelist{this.params.PhaseSimulation.inFileid};
          pathlist = pathlist{this.params.PhaseSimulation.inFileid};
        end
      
        % Initialize Activity:
        act = load(fullfile(pathlist,filelist));
        if iscell(act.act)
          act = act.act{this.params.PhaseSimulation.inCellid};
        else
          act = act.act;
        end
      else
        act = this.params.PhaseSimulation.inActFolder;
      end
      
      if this.params.PhaseSimulation.sampleActOnce
        act = double(rand(size(act)) < act);
      end
      
      if this.numJobs > 1
        this.tempsavepath = fullfile(this.temppath,num2str(this.currJobid));
        mkdir(this.tempsavepath);
      else
        this.tempsavepath = this.temppath;
      end
      
      % Initialize Connections:
      if ischar(this.params.PhaseSimulation.inConnFilename)
        W = load(fullfile(this.workpath, this.params.PhaseSimulation.inConnFilename));
        W = W.W;
      else
        W = this.params.PhaseSimulation.inConnFilename;
      end
      
      if this.params.PhaseSimulation.weightAll ~= 1
        W.w = W.w * this.params.PhaseSimulation.weightAll;
      end
      if this.params.PhaseSimulation.weightExc ~= 1
        W.w(W.w > 0) = W.w(W.w > 0) * this.params.PhaseSimulation.weightExc;
      end
      if this.params.PhaseSimulation.weightInh ~= 1
        if this.params.PhaseSimulation.weightInh==0
          tmpids = W.w < 0;
          W.dx(tmpids) = [];
          W.dy(tmpids) = [];
          W.f0(tmpids) = [];
          W.f1(tmpids) = [];
          W.w(tmpids) = [];
        else
          W.w(W.w < 0) = W.w(W.w < 0) * this.params.PhaseSimulation.weightInh;
        end
      end
      
      % Initialize Phase:
      if ~isempty(this.params.PhaseSimulation.inPhaseFolder)
        inputfolderPhase = fullfile(this.workpath,this.params.PhaseSimulation.inPhaseFolder);
        [ pathlistPhase, filelistPhase ] = dirrec( inputfolderPhase,this.params.PhaseSimulation.inPhaseFilenames );
        if ~isempty(this.params.PhaseSimulation.inFileid)
          filelistPhase = filelistPhase{this.params.PhaseSimulation.inFileid};
          pathlistPhase = pathlistPhase{this.params.PhaseSimulation.inFileid};
        end
        phase = load(fullfile(pathlistPhase,filelistPhase));
        if iscell(phase.phase)
          phase = phase.phase{this.params.PhaseSimulation.inCellid};
        else
          phase = phase.phase;
        end
      elseif ~isempty(this.params.PhaseSimulation.inPhaseFilename)
        phase = load(fullfile(this.workpath,this.params.PhaseSimulation.inPhaseFilename));
        phase = phase.phase;
        phase = phase(:);
      else
        phase = 2*pi*rand(size(act));
        phase = phase(:);
      end
      
      if this.numJobs==1 || (length(this.variableParams)==1 && strcmp(this.variableParams{1}{1},'PhaseSimulation') && strcmp(this.variableParams{1}{2},'inFileid'))
        this.savepath = fullfile(this.workpath,this.params.PhaseSimulation.outPhaseFolder);
      else
        this.savepath = fullfile(this.workpath,this.params.PhaseSimulation.outPhaseFolder, num2str(this.currJobid));
      end
      if exist('inputfolder','var') && exist('pathlist','var')
        this.savepath = fullfile(this.savepath , pathlist(length(inputfolder)+2:end) );
      end
      mkdir(this.savepath);
      
      T = 0:this.params.PhaseSimulation.dt:this.params.PhaseSimulation.tmax;
      fhandle = @(t,x) phaseEvalSparse(this, t, x, act, W);
      fhandleIter = @(t,x,flag) iterFcnSave( this, t, x, act, flag );
      options = odeset('RelTol',10,'AbsTol',0.05,'OutputFcn',fhandleIter ); %RelTol is set to high because phase is anyway circular
      
      %% Start Simulation:
      if strcmp(this.params.PhaseSimulation.odeSolver,'ode1')
        phase = ode1new(fhandle,T,phase(:), options);
      elseif strcmp(this.params.PhaseSimulation.odeSolver,'ode2')
        phase = ode2new(fhandle,T,phase(:), options);
      elseif strcmp(this.params.PhaseSimulation.odeSolver,'ode3')
        phase = ode3new(fhandle,T,phase(:), options);
      elseif strcmp(this.params.PhaseSimulation.odeSolver,'ode4')
        phase = ode4new(fhandle,T,phase(:), options);
      elseif strcmp(this.params.PhaseSimulation.odeSolver,'ode23')
        [~,phase] = ode23(fhandle,T,phase(:), options);
      end
      
      phase = reshape(phase,size(act)); %#ok<NASGU>
      save(fullfile(this.savepath,'phase'),'-v7.3','phase');
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function [ dphase ] = phaseEvalSparse( this, t, phase, activity, W)
      %PHASEEVAL Summary of this function goes here
      %   W is a structure array with elements:
      %     dx,dy,f0,f1,w,delay
      global currNoise;
      
      if this.params.PhaseSimulation.sampleActRepeated
        activity = double(rand(size(activity)) < activity);
      end
      
      disp(num2str(t));
      tmpPhaseSize=size(phase);
      phase = reshape(phase,size(activity)); % x,y,f
      
      numThreads = 1;
      cutsW = round(linspace(1,length(W.w)+1,numThreads+1));
      cutsWend = cutsW(2:end)-1;
      dphasePar = cell(numThreads,1);
      
      for th=1:numThreads
        dphasePar{th} = zeros(size(activity)); % x,y,f
        
        for i=cutsW(th):cutsWend(th)
          %     disp(['i =' num2str(i) ' / ' num2str(cutsWend(th))]);
          
          %% Shift images according to dx:
          if W.dx(i)>=0
            xCutStartIdf0 = W.dx(i);
            xCutStartIdf1 = 0;
            xCutEndIdf0 = 0;
            xCutEndIdf1 = W.dx(i);
          else
            xCutStartIdf0 = 0;
            xCutStartIdf1 = -W.dx(i);
            xCutEndIdf0 = -W.dx(i);
            xCutEndIdf1 = 0;
          end
          
          %% Shift images according to dy:
          if W.dy(i)>=0
            yCutStartIdf0 = W.dy(i);
            yCutStartIdf1 = 0;
            yCutEndIdf0 = 0;
            yCutEndIdf1 = W.dy(i);
          else
            yCutStartIdf0 = 0;
            yCutStartIdf1 = -W.dy(i);
            yCutEndIdf0 = -W.dy(i);
            yCutEndIdf1 = 0;
          end
          
          %% calculate interaction from neuron(x,y,f1) to neuron(x+dx,y+dy,f0)
          phase0 = phase(1+xCutStartIdf0:end-xCutEndIdf0,1+yCutStartIdf0:end-yCutEndIdf0,W.f0(i));
          phase1 = phase(1+xCutStartIdf1:end-xCutEndIdf1,1+yCutStartIdf1:end-yCutEndIdf1,W.f1(i));
          activity0 = activity(1+xCutStartIdf0:end-xCutEndIdf0,1+yCutStartIdf0:end-yCutEndIdf0,W.f0(i));
          activity1 = activity(1+xCutStartIdf1:end-xCutEndIdf1,1+yCutStartIdf1:end-yCutEndIdf1,W.f1(i));
          dphaseTmp = W.w(i) * activity0 .* activity1 .* sin( phase0 - phase1 + this.params.PhaseSimulation.fixedPhaseDelay );
          
          %% phase update:
          dphasePar{th}(1+xCutStartIdf0:end-xCutEndIdf0,1+yCutStartIdf0:end-yCutEndIdf0,W.f0(i)) = ...
            dphasePar{th}(1+xCutStartIdf0:end-xCutEndIdf0,1+yCutStartIdf0:end-yCutEndIdf0,W.f0(i)) - ...
            dphaseTmp;
          
        end
      end
      dphase = sum(cat(4,dphasePar{:}),4);
      clear dphasePar;
      
      %% add global inhibition:
      if ~isempty(this.params.PhaseSimulation.weightInhRadialFcn) && this.params.PhaseSimulation.weightGlobInh > 0
        sumActSinPhase = sum(activity.*sin(phase),3);
        sumActCosPhase = sum(activity.*cos(phase),3);
        
        kernelDim1 = size(activity,1)*2-1;
        kernelDim2 = size(activity,2)*2-1;
        distMat = sqrt(bsxfun(@plus,((1:kernelDim1)-(kernelDim1+1)/2).^2', ((1:kernelDim2)-(kernelDim2+1)/2).^2));
        distMat = feval(this.params.PhaseSimulation.weightInhRadialFcn,distMat);
        distMat = distMat / sum(distMat(:));
        
        sumActSinPhase=conv2(sumActSinPhase,distMat,'same');
        sumActCosPhase=conv2(sumActCosPhase,distMat,'same');
        
        dphaseGlobInh = this.params.PhaseSimulation.weightGlobInh * activity.* ( bsxfun(@times,sin(phase),sumActCosPhase) - bsxfun(@times,cos(phase),sumActSinPhase) );
        
        disp(['mean(abs(dphase)) = ' num2str(mean(abs(dphase(:))))]);
        disp(['mean(abs(dphaseGlobInh)) = ' num2str(mean(abs(dphaseGlobInh(:))))]);
        
        dphase = dphase + dphaseGlobInh;
      end
      
      %% add noise
      if this.params.PhaseSimulation.noiseLevel
        dphase = dphase + this.params.PhaseSimulation.noiseLevel*currNoise;
      end
      
      %% limit the maximum change in phase per iteration:
      maxdphase = this.params.PhaseSimulation.maxdphase;
      numNeuronsOutOfBounds = sum(abs(dphase(:))>maxdphase);
      if numNeuronsOutOfBounds>0
        disp(['numNeuronsOutOfBounds = ' num2str(numNeuronsOutOfBounds) ' = ' num2str(100*numNeuronsOutOfBounds/numel(dphase)) ' %'])
        dphase(dphase>maxdphase) = maxdphase;
        dphase(dphase<-maxdphase) = -maxdphase;
      end
      
      dphase = reshape(dphase,tmpPhaseSize);
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
    
    
    function [ stop ] = iterFcnSave( this, t, phase, activity, flag )
      global currNoise;
            
      if strcmp(flag,'init')
        t = t(1);
      end
      if strcmp(flag,'done')
        return;
      end
      
      phase = reshape(phase,size(activity));

      if this.params.PhaseSimulation.noiseLevel
        if this.params.PhaseSimulation.noiseEMAconst && ~isempty(currNoise)
          currNoise = this.params.PhaseSimulation.noiseEMAconst * currNoise + randn(size(activity));
        else
          currNoise = randn(size(activity));
        end
      end
      
      if mod(t, this.params.PhaseSimulation.saveintervalPhase)==0
        save(fullfile(this.savepath,['phaseIter' num2str(t) '.mat']),'phase');
      end
      
      if mod(t, this.params.PhaseSimulation.saveintervalMeanPhase)==0
        meanPhase = squeeze(circ_mean(phase,ones(size(activity)),3)); %#ok<NASGU>
        save(fullfile(this.savepath,['meanPhaseIter' num2str(t) '.mat']),'meanPhase');
      end
      
      if mod(t, this.params.PhaseSimulation.saveintervalMeanWeightedPhase)==0
        meanWeightedPhase = squeeze(circ_mean(phase,activity,3)); %#ok<NASGU>
        save(fullfile(this.savepath,['meanWeightedPhaseIter' num2str(t) '.mat']),'meanWeightedPhase');
      end
      
      if mod(t, this.params.PhaseSimulation.saveintervalMeanWeightedPhasePng)==0
        meanWeightedPhase = squeeze(circ_mean(phase,activity,3));
        cmap = colormap('hsv');
        meanWeightedPhase = (meanWeightedPhase+pi)/(2*pi);
        meanWeightedPhase = grs2rgb(meanWeightedPhase,cmap,true);
        imwrite(meanWeightedPhase,fullfile(this.savepath,['pngMeanWeightedPhaseIter' num2str(t) '.png']),'png')
      end
      
      if this.params.PhaseSimulation.plotPhase
        sfigure(1);
        imagesc(squeeze(circ_mean(phase,activity,3)));
        colormap hsv;
        caxis([-pi pi]);
        colorbar;
        drawnow;
      end
      
      stop=false;
    end
    
    
  end
  
end

