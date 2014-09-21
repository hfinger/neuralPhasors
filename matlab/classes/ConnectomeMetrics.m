classdef ConnectomeMetrics < Gridjob
  %ConnectomeMetrics Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ConnectomeMetrics(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ConnectomeMetrics.sparse = 0.6;                           % [0.6, ...., 0.75[, graph should not decompose
      this.params.ConnectomeMetrics.heuristics = 1;                         % 1:1:7;
      this.params.ConnectomeMetrics.hScale = 1;                             % vary over scaling of homotopic connections
      
      this.params.ConnectomeMetrics.graphH0 = 1000;                         % vd Heuvel & Sporns (2011) used 1000 permutation networks
      this.params.ConnectomeMetrics.outFilenames = 'ConnectomeMetrics';
      
      % concerning this.params.ConnectomeMetrics.sparse:
        % tic; SPRS = metricsGlobal_wu(SC06, nSamples); toc;                  graph permutation is very expensive
        % Elapsed time is 273.425889 seconds.                                 for densely connected networks
        % tic; FULL = metricsGlobal_wu(SC, nSamples); toc;                    0.7 vs 0.0 sparseness
        % Elapsed time is 1857.124529 seconds.
      
      
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
      
      param = this.params.ConnectomeMetrics;
      
      paths = dataPaths();
      load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat')); % get SC matrix, use average SC
      
      % recalculate avg_ci, only take connections present in >= 75% of subjects to remove outliers?
      % our matrices have way more connections than data used in van den Heuvel & Sporns (2011)
      % resulting in differences of metrics. --> threshold connections/
      
      SC = avg_ci;
      SC(isnan(SC)) = 0;
      SC = SC + SC';                                                        % make SC symmetric --
      % symm.con.: ROI normalization
      % non-symm.con.: row normalization, sparsification
      % trigIds = find(tril(ones([66 66]),0));
 
     %% 
      path = strcat(paths.localTempDir,'/Results/','sp',num2str(param.sparse));
      if ~exist(path,'dir')
          mkdir(path);
      end
      
      SCNorm = normGraph(SC, avg_roi_size, 'ROIprd', false, param.sparse);
      SCMetr = graphMetrics(SCNorm, param.graphH0, path);                                 
      
     %%
      path = strcat(paths.localTempDir,'/Results/','sp',num2str(param.sparse),'/heur',num2str(param.heuristics),'/h',num2str(param.hScale));
      if ~exist(path,'dir')
          mkdir(path);
      end
      
      hSC = homotopHeur(SCNorm, SCMetr, param.heuristics, param.hScale, param.graphH0);
      hSC = normGraph(hSC, avg_roi_size, 'ROIprd', false, param.sparse+42); % sparse > 1 disables modification of sparseness
      hMetr = graphMetrics(hSC, param.graphH0, path);
      hSC = bsxfun(@rdivide,hSC,sum(hSC,2))';                               % norm rows, i.e. input, to sum(CIJ) == 1 

      %%% 
      path = fullfile(this.workpath,param.outFilenames);
      if ~exist(path,'dir')
          mkdir(path);
      end
      savefilename = fullfile(path,num2str(this.currJobid));
      %save([savefilename 'SC.mat' ],'SCNorm','SCMetr');                    % removed: take case hScale = 0 as SC matrix, only when add == true in homotopHeur!     
      save([savefilename 'SC.mat'],'hSC','hMetr');                          % hSC is row-normalized
                  
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
      % loop over 1:this.numJobs to merge mats into 1 struct
      % ...
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

