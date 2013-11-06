classdef Evolution < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

    genId
    
    evoParams
    evoVariableParams
    evoParamComb
    
    individualsAlive = []
    parentsOfIndividuals = [];
    generationOfIndividuals = [];
    
    fitnessAllGen = []
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = Evolution(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.Evolution.paramsAll = struct();
      this.params.Evolution.paramsAllRange = struct(); % specified per structure element as [initMinVal initMaxVal doScaleLog limitRange doWrap]
      this.params.Evolution.numGenerations = 10;
      this.params.Evolution.populationSize = 10;
      this.params.Evolution.numNewIndividualsPerGen = 10;
      this.params.Evolution.initFromOneIndividual = false;
      this.params.Evolution.mutationInitial = 1;
      this.params.Evolution.mutationHalfIter = 3;
      this.params.Evolution.fitnessFcn = []; % has to be specified as a function handle to a function(evoIndividualPath) which returns the fitness 
      this.params.Evolution.startGenId = 1;
      this.params.Evolution.doReplaceLeastFitParent = true;
      
      %%
      this.params.Evolution.startedByIndividual = []; %only for internal use

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
      
      p = this.params.Evolution;
      
      lastGenId = p.startGenId-1;
      lastEvoPath = ['evo' num2str(lastGenId)];
      
      %% calculate mutation rate:
      if ~isempty(p.mutationHalfIter)
        mutationRate = p.mutationInitial / (1 + p.startGenId/p.mutationHalfIter);
      else
        mutationRate = 1/p.startGenId;
      end
      disp(['current mutation rate = ' num2str(mutationRate)])
      
      %% check if there are evolutions of the last generation in the workdir:
      if exist(fullfile(this.workpath,lastEvoPath),'dir')
        
        %% set individual of last generation, which started this job, to finished state:
        if ~isempty(p.startedByIndividual)
          if exist(fullfile(this.workpath,lastEvoPath,'isIndividualFinished'),'dir')
            fclose(fopen(fullfile(this.workpath,lastEvoPath,'isIndividualFinished',num2str(p.startedByIndividual)), 'w'));
          end
        end
        
        %% load old progress
        tmp = load(fullfile(this.temppath,'evoProgress.mat'));
        this.genId = tmp.this.genId;
        this.evoParams = tmp.this.evoParams;
        this.evoVariableParams = tmp.this.evoVariableParams;
        this.evoParamComb = tmp.this.evoParamComb;
        this.individualsAlive = tmp.this.individualsAlive;
        this.parentsOfIndividuals = tmp.this.parentsOfIndividuals;
        this.generationOfIndividuals = tmp.this.generationOfIndividuals;
        this.fitnessAllGen = tmp.this.fitnessAllGen;
        
        %% check if all individuals from last generation are finished:
        numIndividualsInLasGen = sum(this.generationOfIndividuals == lastGenId);
        if exist(fullfile(this.workpath,lastEvoPath,'isIndividualFinished'),'dir')
          filelist = dir(fullfile(this.workpath,lastEvoPath,'isIndividualFinished'));
          filelist = cellfun(@str2num,{filelist.name},'UniformOutput',false);
          filelist = cell2mat(filelist);
          if length(filelist) >= numIndividualsInLasGen
            disp('all individuals of last generation are finished...')
            %try to finish:
            status = rmdir(fullfile(this.workpath,lastEvoPath,'isIndividualFinished'),'s');
            if ~status
              disp('was not able to remove directory isIndividualFinished. Therefore cannot start next generation...')
              return;
            end
          else
            disp('not all individuals of last generation are finished...')
            return;
          end
        else
          disp('there was no isIndividualFinished folder in last generation')
        end
        
        %% compute the fitness of all individuals of the recent simulations
        newIndividualIds = find(this.generationOfIndividuals==lastGenId);
        for evoIndivId=1:length(newIndividualIds)
          individualPath = fullfile(this.workpath,lastEvoPath,['var' num2str(evoIndivId)]);
          if ~isempty(p.fitnessFcn)
            fitness(evoIndivId) = max(0,feval(p.fitnessFcn,individualPath));
            
            %% replace parents:
            if this.parentsOfIndividuals(1,newIndividualIds(evoIndivId)) ~= 0
              fitnessParent1 = this.fitnessAllGen(this.parentsOfIndividuals(1,newIndividualIds(evoIndivId)));
              fitnessParent2 = this.fitnessAllGen(this.parentsOfIndividuals(2,newIndividualIds(evoIndivId)));
              [~,I] = min([fitness(evoIndivId) fitnessParent1 fitnessParent2]);
              if I==1
                this.individualsAlive(newIndividualIds(evoIndivId)) = false;
              elseif I==2
                this.individualsAlive(this.parentsOfIndividuals(1,newIndividualIds(evoIndivId))) = false;
              else
                this.individualsAlive(this.parentsOfIndividuals(2,newIndividualIds(evoIndivId))) = false;
              end
            end
            
          end
        end
        this.fitnessAllGen(newIndividualIds) = fitness;
        
        %% kill individuals to keep population size constant
        numIndividualsToKill = sum(this.individualsAlive) - p.populationSize;
        idsOfAlive = find(this.individualsAlive);
        fitnessOfAlive = this.fitnessAllGen(this.individualsAlive);
        [~,IX] = sort(fitnessOfAlive);
        killIds = idsOfAlive(IX(1:numIndividualsToKill));
        this.individualsAlive(killIds) = false;
        
        %% create new individuals:
        numIndividuals = length(this.fitnessAllGen);
        newIndividuals = [];
        numNewIndividuals = p.numNewIndividualsPerGen;  % p.populationSize - sum(this.individualsAlive);
        idsOfAlive = find(this.individualsAlive);
        fitnessOfAlive = this.fitnessAllGen(idsOfAlive);
        for evoIndivId=1:numNewIndividuals
          %% select from all individuals in all past generations with probability ~ fitness:
          cumsumFitness = [0 cumsum(fitnessOfAlive)];
          parentId1 = idsOfAlive(find( cumsumFitness < rand(1)*cumsumFitness(end) , 1, 'last'));
          parentId2 = idsOfAlive(find( cumsumFitness < rand(1)*cumsumFitness(end) , 1, 'last'));
          weightingOfParent1 = this.fitnessAllGen(parentId1) / (this.fitnessAllGen(parentId1) + this.fitnessAllGen(parentId2));
          
          this.parentsOfIndividuals(:,end+1) = [parentId1, parentId2];
          this.individualsAlive(end+1) = true;
          this.generationOfIndividuals(end+1) = p.startGenId;
        
          individualId = numIndividuals + evoIndivId;
          newIndividuals(end+1) = individualId;
          for job=1:length(p.paramsAllRange)
            this.evoParams(job).params = Evolution.crossoverParams(p.paramsAllRange(job).params,this.evoParams(job).params,parentId1,parentId2,weightingOfParent1);
            this.evoParams(job).params = Evolution.mutateParams(p.paramsAllRange(job).params,this.evoParams(job).params,individualId,mutationRate);
          end
        end
        
      else
        
        %% initialize evolutionary algorithm with random values:
        for job=1:length(p.paramsAllRange)
          if p.initFromOneIndividual
            this.evoParams(job).params = Evolution.initIndividuals(p.paramsAllRange(job).params,p.paramsAll(job).params,p.populationSize);
            for indivId=1:p.populationSize
              this.evoParams(job).params = Evolution.mutateParams(p.paramsAllRange(job).params,this.evoParams(job).params,indivId,mutationRate);
            end
          else
            this.evoParams(job).params = Evolution.selectRandomParams(p.paramsAllRange(job).params,p.populationSize);
          end
        end
        newIndividuals = 1:p.populationSize;
        this.individualsAlive = true(1,p.populationSize);
        this.parentsOfIndividuals = zeros(2,p.populationSize);
        this.generationOfIndividuals = p.startGenId*ones(1,p.populationSize);
        
      end
      
      this.genId = p.startGenId;
      
      %% calculate varParams of evolutionary variables:
      for job=1:length(p.paramsAllRange)
        [this.evoVariableParams{job}, this.evoParamComb{job}] = getVariableParams(this.evoParams(job).params,true);
        
      end
      
      %% update save evolution progress:
      save(fullfile(this.temppath,'evoProgress.mat'),'this');
      
      %% Now start subjobs in subWorkdirs:
      evoPath = ['evo' num2str(this.genId)];
      mkdir(fullfile(this.workpath,evoPath,'isIndividualFinished'));
      for evoIndivId=1:length(newIndividuals)
        
        paramsAll = cell(0);
        
        for job=1:length(p.paramsAll)
          
          if ~isempty(this.evoParamComb{job})
            paramsAll{job} = applyVariableParams(p.paramsAll(job).params, this.evoVariableParams{job}, this.evoParamComb{job}(:,newIndividuals), evoIndivId);
          else
            paramsAll{job} = p.paramsAll(job).params;
          end
          paramsAll{job}.Gridjob.relativeWorkpath = '';
          
          individualPath = fullfile(evoPath,['var' num2str(evoIndivId)]);
          
          %% set workpath and resultpath:
          if this.params.Gridjob.useAbsolutePaths
            paramsAll{job}.Gridjob.useAbsolutePaths = true;
            paramsAll{job}.Gridjob.workpath = fullfile(this.workpath,individualPath);
            paramsAll{job}.Gridjob.resultpath = fullfile(this.resultpath,individualPath);
          elseif ~isempty(this.params.Gridjob.relativeWorkpath)
            paramsAll{job}.Gridjob.relativeWorkpath = fullfile(p.relativeWorkpath,individualPath);
          else
            if iscell(this.constructedFromFolder)
              relativePath=fullfile(this.constructedFromFolder{:});
              paramsAll{job}.Gridjob.relativeWorkpath = fullfile(relativePath,individualPath);
            else
              paramsAll{job}.Gridjob.useAbsolutePaths = true;
              paramsAll{job}.Gridjob.workpath = fullfile(this.constructedFromFolder,individualPath,paramsAll{job}.Gridjob.workpath);
              paramsAll{job}.Gridjob.resultpath = fullfile(this.constructedFromFolder,individualPath,paramsAll{job}.Gridjob.resultpath);
            end
          end
          
        end
        
        %% append Evolution job for next evolution step:
        
        paramsAll{end+1} = this.params;
        paramsAll{end}.Evolution.startGenId = this.genId + 1;
        paramsAll{end}.Evolution.startedByIndividual = evoIndivId;
        paramsAll{end}.Gridjob.initRandStreamWithSeed = rand(1);
        paramsAll{end}.Gridjob.continue = true;
        paramsAll{end}.Gridjob.runLocal = true;
        
        %% start this individual of the next generation:
        disp(['now start individual ' individualPath])
        
        nextJob = Gridjob.createSubclass(paramsAll{1});
        nextJob.paramQueue = paramsAll(2:end);
        nextJob.constructedFromFolder = this.constructedFromFolder;
        start(nextJob);
        clear nextJob;
        
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
  
  methods (Static)
    
    function valStruct = selectRandomParams(rangeStruct,populationSize)
      
      fnames = fieldnames(rangeStruct);
      
      valStruct = struct();
      
      for f=1:length(fnames)
        
        if isstruct(rangeStruct.(fnames{f}))
          
          valStruct.(fnames{f}) = Evolution.selectRandomParams(rangeStruct.(fnames{f}),populationSize);          
          
        else
          
          initMinVal = rangeStruct.(fnames{f})(1);
          initMaxVal = rangeStruct.(fnames{f})(2);
          doScaleLog = rangeStruct.(fnames{f})(3);
          limitRange = rangeStruct.(fnames{f})(4);
          
          if doScaleLog
            
            %% scale values in log scale:
            valuesLogScale = log(initMinVal)+(log(initMaxVal)-log(initMinVal))*rand(1,populationSize);
            values = exp(valuesLogScale);
            valStruct.(fnames{f}) = num2cell(values);
            
          else
            
            %% scale values in linear scale:
            values = initMinVal+(initMaxVal-initMinVal)*rand(1,populationSize);
            valStruct.(fnames{f}) = num2cell(values);
            
          end
          
        end
      
      end
      
    end
    
    function valStruct = initIndividuals(rangeStruct,paramsStruct,populationSize)
      fnames = fieldnames(rangeStruct);
      valStruct = struct();
      for f=1:length(fnames)
        if isstruct(rangeStruct.(fnames{f}))
          valStruct.(fnames{f}) = Evolution.initIndividuals(rangeStruct.(fnames{f}),paramsStruct.(fnames{f}),populationSize);          
        else
          valStruct.(fnames{f}) = repmat( num2cell(paramsStruct.(fnames{f})), 1, populationSize );
        end
      end
    end
    
    function valStruct = mutateParams(rangeStruct,valStruct,individualId,mutationRate)
      
      fnames = fieldnames(rangeStruct);
      
      for f=1:length(fnames)
        
        if isstruct(rangeStruct.(fnames{f}))
          
          valStruct.(fnames{f}) = Evolution.mutateParams(rangeStruct.(fnames{f}),valStruct.(fnames{f}),individualId,mutationRate);          
          
        else
          
          initMinVal = rangeStruct.(fnames{f})(1);
          initMaxVal = rangeStruct.(fnames{f})(2);
          doScaleLog = rangeStruct.(fnames{f})(3);
          limitRange = rangeStruct.(fnames{f})(4);
          doWrap = rangeStruct.(fnames{f})(5);
          
          origVal = valStruct.(fnames{f}){individualId};
          
          if doScaleLog
            origVal = log(origVal);
            initMinVal = log(initMinVal);
            initMaxVal = log(initMaxVal);
          end
          
          range = initMaxVal-initMinVal;
          
          newVal = origVal + range * mutationRate * randn(1);
          
          if limitRange
            if doWrap
              if newVal > initMaxVal
                remaining = newVal - initMaxVal;
                remaining = mod(remaining, range );
                newVal = initMinVal + remaining;
              end
              if newVal < initMinVal
                remaining = initMinVal - newVal;
                remaining = mod(remaining, range );
                newVal = initMaxVal - remaining;
              end
            else
              newVal = min(initMaxVal,newVal);
              newVal = max(initMinVal,newVal);
            end
            assert( newVal>=initMinVal && newVal<=initMaxVal )
          end
          
          if doScaleLog
            newVal = exp(newVal);
          end
          
          valStruct.(fnames{f}){individualId} = newVal;
          
        end
      
      end
      
    end
    
    function valStruct = crossoverParams(rangeStruct,valStruct,parentId1,parentId2,weightingOfParent1)
      
      fnames = fieldnames(rangeStruct);
      
      for f=1:length(fnames)
        
        if isstruct(rangeStruct.(fnames{f}))
          
          valStruct.(fnames{f}) = Evolution.crossoverParams(rangeStruct.(fnames{f}),valStruct.(fnames{f}),parentId1,parentId2,weightingOfParent1);          
          
        else
          
          doScaleLog = rangeStruct.(fnames{f})(3);
          
          origVal1 = valStruct.(fnames{f}){parentId1};
          origVal2 = valStruct.(fnames{f}){parentId2};
          
          if doScaleLog
            origVal1 = log(origVal1);
            origVal2 = log(origVal2);
          end
          
          newVal = origVal1*weightingOfParent1 + origVal2*(1-weightingOfParent1);
          
          if doScaleLog
            newVal = exp(newVal);
          end
          
          valStruct.(fnames{f}){end+1} = newVal;
          
        end
      
      end
      
    end
    
  end
  
end

