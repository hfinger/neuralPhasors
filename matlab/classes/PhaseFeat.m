classdef PhaseFeat < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = PhaseFeat(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.PhaseFeat.inActFolder = 'inAct'; %relative to the workpath
      this.params.PhaseFeat.inActFilenames = 'act.*.mat';
      this.params.PhaseFeat.inPhaseFolder = []; %relative to the workpath
      this.params.PhaseFeat.inPhaseFilenames = [];
      this.params.PhaseFeat.fileid = [];
      this.params.PhaseFeat.outPhaseFeatFolder = 'outAct'; %relative to the workpath
      this.params.PhaseFeat.filtRadius = 2;
      this.params.PhaseFeat.coherenceFactor = 1;
      
      
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
      inputfolder = fullfile(this.workpath,this.params.PhaseFeat.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.PhaseFeat.inActFilenames );
      
      inputfolderPhase = fullfile(this.workpath,this.params.PhaseFeat.inPhaseFolder);
      [ pathlistPhase, filelistPhase ] = dirrec( inputfolderPhase,this.params.PhaseFeat.inPhaseFilenames );
      
      if this.numJobs==1 || (length(this.variableParams)==1 && strcmp(this.variableParams{1}{1},'PhaseFeat') && strcmp(this.variableParams{1}{2},'inFileid'))
        savedir = fullfile(this.workpath,this.params.PhaseFeat.outPhaseFeatFolder);
      else
        savedir = fullfile(this.workpath,this.params.PhaseFeat.outPhaseFeatFolder, num2str(this.currJobid));
      end
      
      if ~isempty(this.params.PhaseFeat.fileid)
        filelist = filelist(this.params.PhaseFeat.fileid);
        pathlist = pathlist(this.params.PhaseFeat.fileid);
        if ~isempty(this.params.PhaseFeat.inPhaseFolder)
          filelistPhase = filelistPhase(this.params.PhaseFeat.fileid);
          pathlistPhase = pathlistPhase(this.params.PhaseFeat.fileid);
        end
      end
      
      for fileid=1:length(filelist)
        disp(fullfile(pathlist{fileid},filelist{fileid}));
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        phase = load(fullfile(pathlistPhase{fileid},filelistPhase{fileid}));
        phase = phase.phase;
        if iscell(act)
          for cellid=1:numel(act)
            feat{cellid} = PhaseFeat.genPhaseFeat(act{cellid}, phase, this.params.PhaseFeat.filtRadius, this.params.PhaseFeat.coherenceFactor);
          end
        else
          feat = PhaseFeat.genPhaseFeat(act, phase, this.params.PhaseFeat.filtRadius, this.params.PhaseFeat.coherenceFactor);
        end
        
        savepath = fullfile(savedir, pathlistPhase{fileid}(length(inputfolderPhase)+2:end));
        mkdir(savepath);
        save(fullfile(savepath,filelist{fileid}),'-v7.3','feat');
        
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
    
    function [ locCoherence ] = genPhaseFeat( activity, phase, filtRadius, coherenceFactor )

      sumLocAct = sum(activity,3);
      if any(sumLocAct(:)==0)
        locAct = bsxfun(@plus,activity,(sumLocAct==0)*1e-10); %add small numbers at points where all activity==0
        sumLocAct = sum(locAct,3);
      end

      sumLocPhase = sum(activity.*exp(coherenceFactor * 1i * phase),3);
      locCoherence = cell(1,length(filtRadius));
      
        locKernel = fspecial('disk',filtRadius);
        sumLocAct = conv2(sumLocAct,locKernel,'same');
        sumLocPhase = conv2(sumLocPhase,locKernel,'same');
        locCoherence = abs(sumLocPhase./sumLocAct);

    end
  end
  
end

