classdef ConnectomeEvaluate < Gridjob
  % ConnectomeEvaluate: class to get plots and statistics
  %                     from post-simulation evaluation
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ConnectomeEvaluate(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ConnectomeEvaluate.sparse = 0.6;                           % [0.6, ...., 0.75[, graph should not decompose
      this.params.ConnectomeEvaluate.heuristics = 1;                         % 1:1:7;
      this.params.ConnectomeEvaluate.hScale = 1;                             % vary over scaling of homotopic connections
      
      this.params.ConnectomeEvaluate.graphH0 = 1000;                         % vd Heuvel & Sporns (2011) used 1000 permutation networks
      this.params.ConnectomeEvaluate.outFilenames = 'ConnectomeMetrics';
      
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
      
      param = this.params.ConnectomeEvaluate;
      
      paths = dataPaths();
      load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat')); % get SC matrix, use average SC
      
      %% descriptive work
      
      % evaluate global and local prediction errors,
      % do h vs k plots
      
      % plots of cortexTopology <-- which matrix?
      
      cortexTopology(SC, SCMetrics,sub)
      
      %%
      
      % gradient descent: get 1 new matrix eSC, (free of parameters h, k)
      % obtain graph metrices of this matrix
      
      
      % get performance matrices, plot h vs k
      SimCmp = load(fullfile(paths.workdir,'/pebel','/20150414_SAR_Metrics','/CompareWithEEG/','compareSimExp5.mat'));
      SimCmp.overFreq.coh.rho % SpxHeur x h x k
      
      for sp = 1:2
        for h = 1:2
          figure();
          pMat = squeeze(SimCmp.overFreq.coh.rho(sp,h,:,:));
          imagesc(pMat);
          [M,I] = max(pMat(:));
          [hi, ki] = ind2sub(size(pMat), I);
          title(strcat('Perf. ', num2str(M) , ' for h=',num2str(hi),' k=',num2str(ki)))
        end
      end
      
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

