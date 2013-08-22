classdef ActPhaseCombine < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ActPhaseCombine(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ActPhaseCombine.inActFolder = 'inAct'; %relative to the workpath
      this.params.ActPhaseCombine.inActFilenames = 'act.*.mat';
      this.params.ActPhaseCombine.inPhaseFolder = []; %relative to the workpath
      this.params.ActPhaseCombine.inPhaseFilenames = [];
      this.params.ActPhaseCombine.fileid = [];
      this.params.ActPhaseCombine.outActFolder = 'outAct'; %relative to the workpath
      
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

      
      savedir = fullfile(this.workpath,this.params.ActPhaseCombine.outActFolder);
      
      if ~isempty(this.params.ActPhaseCombine.fileid)
        filelist = filelist(this.params.ActPhaseCombine.fileid);
        pathlist = pathlist(this.params.ActPhaseCombine.fileid);
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
              act{cellid} = act{cellid} .* exp(1i * phase{cellid});
            end
          end
        else
          if ~isempty(this.params.ActPhaseCombine.inPhaseFolder)
            act = act .* exp(1i * phase);
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
  
end

