classdef ActPhaseSplit < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ActPhaseSplit(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ActPhaseSplit.inActPhaseFolder = 'inAct'; %relative to the workpath
      this.params.ActPhaseSplit.inActPhaseFilenames = 'act.*.mat';
      this.params.ActPhaseSplit.fileid = [];
      this.params.ActPhaseSplit.outActFolder = 'outAct'; %relative to the workpath
      this.params.ActPhaseSplit.outPhaseFolder = 'outPhase'; %relative to the workpath
      
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
      inputfolder = fullfile(this.workpath,this.params.ActPhaseSplit.inActPhaseFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.ActPhaseSplit.inActPhaseFilenames );

      savedirAct = fullfile(this.workpath,this.params.ActPhaseSplit.outActFolder);
      savedirPhase = fullfile(this.workpath,this.params.ActPhaseSplit.outPhaseFolder);
      
      if ~isempty(this.params.ActPhaseSplit.fileid)
        filelist = filelist(this.params.ActPhaseSplit.fileid);
        pathlist = pathlist(this.params.ActPhaseSplit.fileid);
      end
      
      for fileid=1:length(filelist)
        disp(fullfile(pathlist{fileid},filelist{fileid}));
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        if iscell(act)
          phase = cell(size(act));
          for cellid=1:numel(act)
            if ~isempty(this.params.ActPhaseSplit.inActPhaseFolder)
              phase{cellid} = angle(act{cellid});
              act{cellid} = abs(act{cellid});
            end
          end
        else
          if ~isempty(this.params.ActPhaseSplit.inActPhaseFolder)
            phase = angle(act);
            act = abs(act);
          end
        end
        if ~isempty(this.params.ActPhaseSplit.outActFolder)
          savepathAct = fullfile(savedirAct, pathlist{fileid}(length(inputfolder)+2:end));
          mkdir(savepathAct);
          save(fullfile(savepathAct,filelist{fileid}),'-v7.3','act');
        end
        if ~isempty(this.params.ActPhaseSplit.outPhaseFolder)
          savepathPhase = fullfile(savedirPhase, pathlist{fileid}(length(inputfolder)+2:end));
          mkdir(savepathPhase);
          save(fullfile(savepathPhase,filelist{fileid}),'-v7.3','phase');
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
  
end

