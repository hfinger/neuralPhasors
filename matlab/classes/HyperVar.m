classdef HyperVar < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
    
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = HyperVar(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.HyperVar.paramsAll = struct();
      this.params.HyperVar.paramsAllHyperVars = struct();
      this.params.HyperVar.subjobWorkdir = 'hypervars';
      
      %% only for internal use:
      this.params.HyperVar.startedByHypervarId = [];
      
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    %% Start: is executed before all individual parameter jobs are started
    function startJob(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
      
      disp('this excecutes before the parallel job is started');
      
      %%%% END EDIT HERE:                                        %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      p = this.params.HyperVar;
      
      hyperVariableParams = cell(1,length(p.paramsAll));
      hyperParamComb = cell(1,length(p.paramsAll));
      for job=1:length(p.paramsAll)
        [hyperVariableParams{job}, hyperParamComb{job}] = getVariableParams(p.paramsAllHyperVars(job).params,true);
      end
      
      numHyperVars = size(hyperParamComb{1},2);
      
      subjobWorkdir = fullfile(this.workpath,p.subjobWorkdir);
      
      if isempty(p.startedByHypervarId)
        
        %% Now start subjobs in subWorkdirs:
        mkdir(fullfile(subjobWorkdir,'isHypervarFinished'));
        
        for hyperVarId=1:length(numHyperVars)
          
          paramsAll = cell(0);
          
          for job=1:length(p.paramsAll)
            
            if ~isempty(hyperParamComb{job})
              paramsAll{job} = applyVariableParams(p.paramsAll(job).params, hyperVariableParams{job}, hyperParamComb{job}, hyperVarId);
            else
              paramsAll{job} = p.paramsAll(job).params;
            end
            paramsAll{job}.Gridjob.relativeWorkpath = '';
            
            hyperVarPath = fullfile(subjobWorkdir,['var' num2str(hyperVarId)]);
            
            %% set workpath and resultpath:
            if this.params.Gridjob.useAbsolutePaths
              paramsAll{job}.Gridjob.useAbsolutePaths = true;
              paramsAll{job}.Gridjob.workpath = fullfile(this.workpath,hyperVarPath);
              paramsAll{job}.Gridjob.resultpath = fullfile(this.resultpath,hyperVarPath);
            elseif ~isempty(this.params.Gridjob.relativeWorkpath)
              paramsAll{job}.Gridjob.relativeWorkpath = fullfile(p.relativeWorkpath,hyperVarPath);
            else
              if iscell(this.constructedFromFolder)
                relativePath=fullfile(this.constructedFromFolder{:});
                paramsAll{job}.Gridjob.relativeWorkpath = fullfile(relativePath,hyperVarPath);
              else
                paramsAll{job}.Gridjob.useAbsolutePaths = true;
                paramsAll{job}.Gridjob.workpath = fullfile(this.constructedFromFolder,hyperVarPath,paramsAll{job}.Gridjob.workpath);
                paramsAll{job}.Gridjob.resultpath = fullfile(this.constructedFromFolder,hyperVarPath,paramsAll{job}.Gridjob.resultpath);
              end
            end
            
          end
          
          %% append HyperVar job for next HyperVar step:
          
          paramsAll{end+1} = this.params;
          paramsAll{end}.HyperVar.startedByHypervarId = hyperVarId;
          paramsAll{end}.Gridjob.initRandStreamWithSeed = rand(1);
          paramsAll{end}.Gridjob.continue = true;
          paramsAll{end}.Gridjob.runLocal = true;
          
          %% start this individual of the next generation:
          disp(['now start individual ' hyperVarPath])
          
          nextJob = Gridjob.createSubclass(paramsAll{1});
          nextJob.paramQueue = paramsAll(2:end);
          nextJob.constructedFromFolder = this.constructedFromFolder;
          
          this.paramQueue = []; %% do not start next jobs directly.
          
          start(nextJob);
          clear nextJob;
          
        end
        
        
        
        
      else
        
        %% set individual of last generation, which started this job, to finished state:
        if exist(fullfile(this.workpath,subjobWorkdir,'isHypervarFinished'),'dir')
          fclose(fopen(fullfile(this.workpath,subjobWorkdir,'isHypervarFinished',num2str(p.startedByHypervarId)), 'w'));
        end
        
        %% check if all HyperVars are finished:
        if exist(fullfile(this.workpath,subjobWorkdir,'isHypervarFinished'),'dir')
          filelist = dir(fullfile(this.workpath,subjobWorkdir,'isHypervarFinished'));
          filelist = cellfun(@str2num,{filelist.name},'UniformOutput',false);
          filelist = cell2mat(filelist);
          if length(filelist) >= numHyperVars
            disp('all hypervars are finished...')
            %try to finish:
            status = rmdir(fullfile(this.workpath,subjobWorkdir,'isHypervarFinished'),'s');
            if ~status
              disp('was not able to remove directory isHypervarFinished. Therefore cannot start next generation...')
              return;
            end
          else
            disp('not all individuals of last generation are finished...')
            return;
          end
        else
          disp('there was no isHypervarFinished folder')
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

