classdef PreprocessZtrafo < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = PreprocessZtrafo(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.PreprocessZtrafo.inActFolder = 'labelMeInput'; %relative to the workpath
      this.params.PreprocessZtrafo.inActFilenames = 'act.*.mat';
      this.params.PreprocessZtrafo.outWeightsFolder = 'zTrafoWeights'; %relative to the workpath
      this.params.PreprocessZtrafo.inNumChannels = 3;
      
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
      inputfolder = fullfile(this.workpath,this.params.PreprocessZtrafo.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.PreprocessZtrafo.inActFilenames );
      
      disp('now compute the mean...')
      tmpSumAct = 0;
      tmpSumActSquare = 0;
      tmpNumAdded = 0;
      
      for fileid=1:numel(filelist)
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        if iscell(act)
          for cellid = 1:length(act)
            tmpSumAct = tmpSumAct + sum(reshape(act{cellid},[size(act{cellid},1)*size(act{cellid},2) size(act{cellid},3)]),1);
            tmpSumActSquare = tmpSumActSquare + sum(reshape(act{cellid}.^2,[size(act{cellid},1)*size(act{cellid},2) size(act{cellid},3)]),1);
            tmpNumAdded = tmpNumAdded + size(act{cellid},1)*size(act{cellid},2);
          end
        else
            tmpSumAct = tmpSumAct + sum(reshape(act,[size(act,1)*size(act,2) size(act,3)]),1);
            tmpSumActSquare = tmpSumActSquare + sum(reshape(act.^2,[size(act,1)*size(act,2) size(act,3)]),1);
            tmpNumAdded = tmpNumAdded + size(act,1)*size(act,2);
        end
      end
      
      meanAct = tmpSumAct / tmpNumAdded;
      meanActSquare = tmpSumActSquare / tmpNumAdded;
      stdAct = sqrt(meanActSquare - meanAct.^2);
      
      conn.inputSubtract = reshape(meanAct, [1 1 length(stdAct(:))]);
      conn.inputScaling = reshape(1./stdAct, [1 1 length(stdAct(:))]);

      savedir = fullfile(this.workpath,this.params.PreprocessZtrafo.outWeightsFolder);
      if this.numJobs > 1
        savedir = fullfile(savedir,num2str(this.currJobid));
      end
      mkdir(savedir);
      
      save(fullfile(savedir,'weights.mat'),'-struct','conn');
      
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

