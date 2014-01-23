classdef Gridjob
  
  properties
    params
    paramQueue
    currJobid
    numJobs
    constructedFromFolder
    workpath
    resultpath
    temppath
    variableParams
    paramComb
    isfinishedDir
  end
  
  methods
    
    %% Constructor:
    function this = Gridjob(varargin)
      
      % set the standard parameters:
      this.params.Gridjob.runLocal = true;
      this.params.Gridjob.restartable = false;
      this.params.Gridjob.continue = false;
      this.params.Gridjob.exclusive = false;
      this.params.Gridjob.combParallel = false;
      this.params.Gridjob.requiremf = []; % in MB
      this.params.Gridjob.wc_host = '!ramsauer';%'ananke|calvin|carpo|daphne|elara|erato|fred|hobbes|isonoe|jupiter|klio|kore|leda|mars|melete|mneme|neptune|saturn|shaggy|thebe|urania|velma|venus'; %|sinope
      this.params.Gridjob.queue = []; % should be empty or 'nbp.q' or 'all.q'
      this.params.Gridjob.jobname = 'job';
      this.params.Gridjob.useAbsolutePaths = false;
      this.params.Gridjob.workpath = '';
      this.params.Gridjob.resultpath = '';
      this.params.Gridjob.deleteTmpFolder = false;
      this.params.Gridjob.initRandStreamWithJobid = false;
      this.params.Gridjob.initRandStreamWithSeed = [];
      this.params.Gridjob.fhandleFinish = [];
      this.params.Gridjob.remoteStart = false;
      this.params.Gridjob.newGrid = true;
      this.params.Gridjob.runOnHPC = false;
      this.params.Gridjob.requiredThreads = '4';%'2-8';
      this.params.Gridjob.matlabpool = 0;
      this.params.Gridjob.relativeWorkpath = [];
      this.params.Gridjob.walltime = 3600;
      
      % save folder from which the object is constructed
      this.constructedFromFolder = pwd;
      paths = dataPaths( );
      %find relative path from pwd to paths.paramdir
      filesepTmp=filesep();
      if filesepTmp=='\'
        filesepTmp='\\';
      end
      partsParamdir = strread(fullfile(paths.paramdir),'%s','delimiter',filesepTmp);
      partsPwd = strread(pwd,'%s','delimiter',filesepTmp);
      matching = true;
      for k=1:min(length(partsParamdir),length(partsPwd))
        if ~strcmp(partsParamdir{k},partsPwd{k})
          matching = false;
          break;
        end
      end
      if matching
        this.constructedFromFolder = partsPwd(length(partsParamdir)+1:end);
      else
        this.constructedFromFolder = pwd;
      end
      
      this = this.init(varargin{:});
    end
    
    
    function this = init(this,initParams)
      % initParams should be either of:
      % 1) a param struct
      % 2) a cell of param structs
      % 3) a string specifying the path to jobDesc.mat
      
      if nargin>1
        if ischar(initParams)
          classProperties = load(initParams);
          if isfield(classProperties,'this')
            fnames = fieldnames(classProperties.this);
            for i=1:length(fnames)
              if isstruct(this.(fnames{i})) && isstruct(classProperties.this.(fnames{i}))
                this.(fnames{i}) = updateStruct(this.(fnames{i}),classProperties.this.(fnames{i}));
              else
                this.(fnames{i}) = classProperties.this.(fnames{i});
              end
            end
          else
            fnames = fieldnames(classProperties);
            for i=1:length(fnames)
              if isfield(this,fnames{i})
                if isstruct(this.(fnames{i})) && isstruct(classProperties.(fnames{i}))
                  this.(fnames{i}) = updateStruct(this.(fnames{i}),classProperties.(fnames{i}));
                else
                  this.(fnames{i}) = classProperties.(fnames{i});
                end
              end
            end
          end
        else
          if iscell(initParams)
            this.paramQueue = initParams(2:end);
            initParams = initParams{1};
          end
          this.params = updateStruct(this.params,initParams);
        end
      end
      
    end
    
    
    %% Start all the jobs:
    function this = start(this)
      
      paths = dataPaths( );
      
      %if this is the gridjob class, then create a subclass with correct type
      if strcmp(class(this),'Gridjob')
        tmp = this.paramQueue;
        this = Gridjob.createSubclass(this.params);
        this.paramQueue = tmp;
%         start(this);
%         return;
      end
      
      %create variableParams and paramComb (to describe parameter sweeps):
      [this.variableParams, this.paramComb] = getVariableParams(this.params,this.params.Gridjob.combParallel);
      this.numJobs = size(this.paramComb,2);
      if this.numJobs == 0
        this.numJobs = 1;
      end
      
      % set workpath and resultpath:
      if this.params.Gridjob.useAbsolutePaths
        localWorkpath = this.params.Gridjob.workpath;
        this.workpath = this.params.Gridjob.workpath;
        this.resultpath = this.params.Gridjob.resultpath;
      elseif ~isempty(this.params.Gridjob.relativeWorkpath)
        relativePath=this.params.Gridjob.relativeWorkpath;
        localWorkpath = fullfile(paths.workdir,relativePath,this.params.Gridjob.workpath);
        if this.params.Gridjob.runLocal
          this.workpath = fullfile(paths.workdir,relativePath,this.params.Gridjob.workpath);
          this.resultpath = fullfile(paths.resultsdir,relativePath,this.params.Gridjob.resultpath);
        else
          this.workpath = fullfile(paths.sge_workdir,relativePath,this.params.Gridjob.workpath);
          this.resultpath = fullfile(paths.sge_resultsdir,relativePath,this.params.Gridjob.resultpath);
        end
      else
        if iscell(this.constructedFromFolder)
          relativePath=fullfile(this.constructedFromFolder{:});
          localWorkpath = fullfile(paths.workdir,relativePath,this.params.Gridjob.workpath);
          if this.params.Gridjob.runLocal
            this.workpath = fullfile(paths.workdir,relativePath,this.params.Gridjob.workpath);
            this.resultpath = fullfile(paths.resultsdir,relativePath,this.params.Gridjob.resultpath);
          else
            this.workpath = fullfile(paths.sge_workdir,relativePath,this.params.Gridjob.workpath);
            this.resultpath = fullfile(paths.sge_resultsdir,relativePath,this.params.Gridjob.resultpath);          
          end
        else
          localWorkpath = fullfile(this.constructedFromFolder,this.params.Gridjob.workpath);
          this.workpath = fullfile(this.constructedFromFolder,this.params.Gridjob.workpath);
          this.resultpath = fullfile(this.constructedFromFolder,this.params.Gridjob.resultpath);
        end
      end
      
      this.temppath = fullfile(this.workpath,['temp_' this.params.Gridjob.jobname]);
      localTemppath = fullfile(localWorkpath,['temp_' this.params.Gridjob.jobname]);

      %create directories:
      if exist(localTemppath, 'dir') && ~this.params.Gridjob.continue
        rmdir(localTemppath,'s');
      end
      mkdir(localTemppath);
      mkdir(fullfile(localTemppath,'isfinished'));
      if ~exist(localWorkpath, 'dir')
        mkdir(localWorkpath);
      end
      
      if sum(strcmp('startJob',methods(this)))
        startJob(this);
      end
      
      if ~this.params.Gridjob.runLocal
        this.temppath = strrep(this.temppath,'\','/');
        this.workpath = strrep(this.workpath,'\','/');
        this.resultpath = strrep(this.resultpath,'\','/');
      end
      
      %save the job description to the temp folder:
      jobDescPath = fullfile(localTemppath,'jobDesc.mat');
      propNames = properties(this);
      for i=1:length(propNames)
        jobDesc.(propNames{i}) = this.(propNames{i}); %#ok<STRNU>
      end
      save(jobDescPath,'-struct','jobDesc');
%       clear jobDesc;
      
      if this.params.Gridjob.runLocal
        
        %start jobs sequentially in this matlab instance
        for jobid=1:this.numJobs
%           Gridjob.startJobid(jobDescPath,jobid);
          Gridjob.startJobid(jobDesc,jobid);
        end
        
      elseif this.params.Gridjob.runOnHPC
        
        if isfield(paths,'hpc_pathToAddScriptPaths')
          pathToAddScriptPaths = paths.hpc_pathToAddScriptPaths;
        else
          pathToAddScriptPaths = which('addScriptPaths');
          pathToAddScriptPaths = pathToAddScriptPaths(1:end-17);
        end
        
        %generate job script:
        jobscriptpath=fullfile(localTemppath,[this.params.Gridjob.jobname '.sh']);
        jobscriptpathRemote=[this.temppath '/' this.params.Gridjob.jobname '.sh'];
        fid=fopen(jobscriptpath, 'w'); 
        fprintf(fid,'#!/bin/bash\n');
        fprintf(fid,'#PBS -J 1-%u\n',this.numJobs);
        fprintf(fid,'#PBS -N %s\n',this.params.Gridjob.jobname);
        fprintf(fid,'#PBS -l ncpus=%s\n',this.params.Gridjob.requiredThreads);
        fprintf(fid,'#PBS -l walltime=%u\n',this.params.Gridjob.walltime);
        fprintf(fid,'#PBS -q workq\n');
        
        fprintf(fid,['/sw/sdev/matlab/R2013b/bin/matlab -nodisplay -r "'...
          'cd ' this.workpath '; addpath(' char(39) pathToAddScriptPaths char(39) '); addScriptPaths(); '...
          'Gridjob.startJobid(' char(39) this.temppath '/jobDesc.mat' char(39) ',$SGE_TASK_ID);" -c ~/ikw_license.dat\n']);
        fprintf(fid,'exit 0');
        fclose(fid);
        
        %check if qsub command is missing
        if system('command -v qsub >/dev/null 2>&1')
          qsubMissing = true;
        else
          qsubMissing = false;
        end
        
        %now start the grid job:
        if qsubMissing || this.params.Gridjob.remoteStart
          if ispc
            systemCmd = ['plink.exe ' paths.plinkArg ' "qsub ' jobscriptpathRemote '"'];
            disp('Execute:');
            disp(systemCmd);
            system(systemCmd);
          else
            system(['ssh isonoe "cd ' pwd '; qsub ' jobscriptpathRemote '"']);
          end
        else
          system(['qsub ' jobscriptpath]);
        end
        
      else
        
        if isfield(paths,'sge_pathToAddScriptPaths')
          pathToAddScriptPaths = paths.sge_pathToAddScriptPaths;
        else
          pathToAddScriptPaths = which('addScriptPaths');
          pathToAddScriptPaths = pathToAddScriptPaths(1:end-17);
        end
        
        %generate job script:
        jobscriptpath=fullfile(localTemppath,[this.params.Gridjob.jobname '.sh']);
        jobscriptpathRemote=[this.temppath '/' this.params.Gridjob.jobname '.sh'];
        fid=fopen(jobscriptpath, 'w'); 
        fprintf(fid,'#!/bin/bash\n');
        fprintf(fid,'#$ -t 1:%u\n',this.numJobs);
        fprintf(fid,'#$ -N %s\n',this.params.Gridjob.jobname);
        if this.params.Gridjob.exclusive
          fprintf(fid,'#$ -l excl=true\n');
        end
        if this.params.Gridjob.newGrid
          if ~isempty(this.params.Gridjob.requiremf)
            fprintf(fid,'#$ -l mem=%uM\n',this.params.Gridjob.requiremf);
          else
            fprintf(fid,'#$ -l mem=%uM\n',5000);
          end
          fprintf(fid,'#$ -pe matlab %s\n',this.params.Gridjob.requiredThreads);
          fprintf(fid,'#$ -m n\n');
          
          if ~isempty(this.params.Gridjob.wc_host)
            fprintf(fid,['#$ -l h=' this.params.Gridjob.wc_host '\n']);
          end
          if ~isempty(this.params.Gridjob.queue)
            fprintf(fid,['#$ -q ' this.params.Gridjob.queue '\n']);
          end
        else
          if ~isempty(this.params.Gridjob.requiremf)
            fprintf(fid,'#$ -l mf=%um\n',this.params.Gridjob.requiremf);
          end
          fprintf(fid,'#$ -l lic_matlab=1\n');

          if isempty(this.params.Gridjob.wc_host)
            fprintf(fid,'#$ -q ikw@*&!(*ramsauer*|*hobbes*|*shaggy*)\n');
          else
            fprintf(fid,['#$ -q ikw@' this.params.Gridjob.wc_host '\n']);
          end
        end
        
        fprintf(fid,'#$ -wd %s\n',this.temppath);
        fprintf(fid,['matlab -nodisplay -r "'...
          'addpath(' char(39) pathToAddScriptPaths char(39) '); addScriptPaths(); '...
          'Gridjob.startJobid(' char(39) this.temppath '/jobDesc.mat' char(39) ',$SGE_TASK_ID);"\n']);
        fprintf(fid,'exit 0');
        fclose(fid);
        
        %check if qsub command is missing
        if system('command -v qsub >/dev/null 2>&1')
          qsubMissing = true;
        else
          qsubMissing = false;
        end
        
        %now start the grid job:
        if qsubMissing || this.params.Gridjob.remoteStart
          if ispc
%             system(['putty.exe -ssh -2 -m c:"cd ' pwd '; qsub ' jobscriptpath '" shaggy']);
%             systemCmd = ['plink.exe -ssh -i ' paths.puttyPrivateKey ' ' paths.sshUsername '@' paths.sshServer ' "qsub ' jobscriptpathRemote '"'];
            systemCmd = ['plink.exe ' paths.plinkArg ' "qsub ' jobscriptpathRemote '"'];
            disp('Execute:');
            disp(systemCmd);
            system(systemCmd);
          else
            system(['ssh isonoe "cd ' pwd '; qsub ' jobscriptpathRemote '"']);
          end
        else
          system(['qsub ' jobscriptpath]);
        end
      end
    end
    
    %% overload this method with the algorithm
    function run(this) %#ok<MANU>
      %...
    end
    
    %% finishJobs and clean up
    function finishJobs(this)
      if this.params.Gridjob.deleteTmpFolder
        rmdir(this.temppath,'s');
      end
      
      if ~isempty(this.params.Gridjob.fhandleFinish)
        tmpPwd = pwd;
        cd(this.workpath);
        feval(this.params.Gridjob.fhandleFinish);
        cd(tmpPwd);
      end
      
      %% start next job:
      if ~isempty(this.paramQueue)
        nextJob = Gridjob.createSubclass(this.paramQueue{1});
        nextJob.paramQueue = this.paramQueue(2:end);
        nextJob.constructedFromFolder = this.constructedFromFolder;
        start(nextJob);
      end
      
    end
    
    
  end
  
  methods (Static)
    
    %% start the individual parameter job
    function startJobid(jobDescPath,jobid)
      if ischar(jobDescPath)
        jobDesc = load(jobDescPath);
      else
        jobDesc = jobDescPath;
      end
      
      this = Gridjob.createSubclass(jobDesc.params);
%       this = updateStruct(this,jobDesc);
      propNames = properties(this);
      for i=1:length(propNames)
        if isfield(jobDesc,propNames{i})
          if isstruct(this.(propNames{i}))
            this.(propNames{i}) = updateStruct(this.(propNames{i}),jobDesc.(propNames{i}));
          else
            this.(propNames{i}) = jobDesc.(propNames{i});
          end
        end
      end
      
      this.currJobid = jobid;
      this.params = applyVariableParams(this.params,this.variableParams,this.paramComb,jobid);
      
%       numSlots = getenv('NSLOTS');
%       numSlots = round(str2double(numSlots));
%       disp(['numSlots reserved: ' num2str(numSlots)]);
%       disp(['numSlots matlab had before: ' num2str(maxNumCompThreads())]);
%       if ~isnan(numSlots)
%         maxNumCompThreads(numSlots);
%       end

      if this.params.Gridjob.matlabpool>0
        matlabpool(this.params.Gridjob.matlabpool);
      end
%       disp(['numSlots matlab has: ' num2str(maxNumCompThreads())]);
      
      if this.params.Gridjob.initRandStreamWithJobid
        if sum(strcmp('setGlobalStream',methods('RandStream')))
          RandStream.setGlobalStream(RandStream('mt19937ar','Seed',jobid));
        else
          RandStream.setDefaultStream(RandStream('mt19937ar','Seed',jobid));
        end
      end
      if ~isempty(this.params.Gridjob.initRandStreamWithSeed)
        if sum(strcmp('setGlobalStream',methods('RandStream')))
          RandStream.setGlobalStream(RandStream('mt19937ar','Seed',this.params.Gridjob.initRandStreamWithSeed));
        else
          RandStream.setDefaultStream(RandStream('mt19937ar','Seed',this.params.Gridjob.initRandStreamWithSeed));
        end
      end
      
      disp(['start jobid: ' num2str(jobid)]);
      run(this);
      fclose(fopen(fullfile(this.temppath,'isfinished',num2str(jobid)), 'w'));
      
      %check if all jobs are finished:
      filelist = dir(fullfile(this.temppath,'isfinished'));
      filelist = cellfun(@str2num,{filelist.name},'UniformOutput',false);
      filelist = cell2mat(filelist);
      if length(filelist) >= this.numJobs
        disp('all jobs are finished')
        %try to finish:
        status = rmdir(fullfile(this.temppath,'isfinished'),'s');
        if status
          finishJobs(this);
        end
      end
    end
    
    function new = createSubclass(params)
      fnames=fieldnames(params);
      for k=1:length(fnames)
        if ~strcmp(fnames{k},'Gridjob')
          break;
        end
      end
      new = feval(fnames{k},params);
    end
    
    
  end
  
   
end