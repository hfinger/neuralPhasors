classdef ActConcat < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ActConcat(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ActConcat.inActFolder1 = 'inAct'; %relative to the workpath
      this.params.ActConcat.inActFilenames1 = 'act.*.mat';
      this.params.ActConcat.inActFolder2 = []; %relative to the workpath
      this.params.ActConcat.inActFilenames2 = [];
      this.params.ActConcat.fileid = [];
      this.params.ActConcat.outActFolder = 'outAct'; %relative to the workpath
      
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
      inputfolder = fullfile(this.workpath,this.params.ActConcat.inActFolder1);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.ActConcat.inActFilenames1 );

      inputfolder2 = fullfile(this.workpath,this.params.ActConcat.inActFolder2);
      [ pathlist2, filelist2 ] = dirrec( inputfolder2,this.params.ActConcat.inActFilenames2 );

      
      savedir = fullfile(this.workpath,this.params.ActConcat.outActFolder);
      
      if ~isempty(this.params.ActConcat.fileid)
        filelist = filelist(this.params.ActConcat.fileid);
        pathlist = pathlist(this.params.ActConcat.fileid);
        filelist2 = filelist2(this.params.ActConcat.fileid);
        pathlist2 = pathlist2(this.params.ActConcat.fileid);
      end
      
      for fileid=1:length(filelist)
        disp(fullfile(pathlist{fileid},filelist{fileid}));
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        act2 = load(fullfile(pathlist2{fileid},filelist2{fileid}));
        act2 = act2.act;
        if iscell(act)
          for cellid=1:numel(act)
            act{cellid} = cat(3, act{cellid} , act2{cellid});
          end
        else
          act = cat(3, act, act2);
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

