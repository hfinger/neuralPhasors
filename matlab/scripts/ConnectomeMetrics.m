classdef ConnectomeMetrics < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
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
      
      this.params.ConnectomeMetrics.graphH0 = true;                         % flag to enable/disable generation of permutation graphs
      this.params.ConnectomeMetrics.outFilenames = 'results';
      
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
      
      load('/net/store/nbp/phasesim/databases/SC_Bastian/dti_20141209_preprocessed.mat');       % get SC matrix, use average SC
      
      % recalculate avg_ci, only take connections present in >= 75% of subjects to remove outliers?
      % our matrices have way more connections than data used in van den Heuvel & Sporns (2011)
      % resulting in differences of metrics. --> threshold connections/
      
      % mat = zeros(66,66);
      % for i = 1:22
      %     if (isempty(ci{1,i})) continue; end
      %     ci{1,i}(isnan(ci{1,i}))=0;
      %     mat = mat + ci{1,i};
      %
      % end
      % mat = mat/(22-3);
      
      
      SC = avg_ci;
      SC(isnan(SC)) = 0;
      SC = SC + SC';                                                                             % make SC symmetric --
      % symm.con.: ROI normalization
      % non-symm.con.: row normalization, sparsification
      % trigIds = find(tril(ones([66 66]),0));
 
     %% 
      path = strcat('Results/sp',num2str(this.params.ConnectomeMetrics.sparse));
      if ~exist(path,'dir')
          mkdir(path);
      end
      
      SCNorm = normGraph(SC, avg_roi_size, 'ROIsum', false, this.params.ConnectomeMetrics.sparse);
      SCMetr = graphMetrics(SCNorm, path);
      
     %%
      path = strcat('Results/sp',num2str(this.params.ConnectomeMetrics.sparse), '/heur',num2str(this.params.ConnectomeMetrics.heuristics), '/h',num2str);
      if ~exist(path,'dir')
          mkdir(path);
      end
      
      hSC = homotopHeur(SCNorm, SCMetr, avg_dist, this.params.ConnectomeMetrics.heuristics, this.params.ConnectomeMetrics.hScale);
      hSC = normGraph(hSC, avg_roi_size, 'ROIsum', false, this.params.ConnectomeMetrics.sparse);
      hMetr = graphMetrics(hSC, path);
      hSC = bsxfun(@rdivide,hSC,sum(hSC,2))';                               % norm rows, i.e. input, to sum(CIJ) == 1 

      
      %%% 
      path = fullfile(this.workpath,this.params.ConnectomeMetrics.outFilenames);
      if ~exist(path,'dir')
          mkdir(path);
      end
      savefilename = fullfile(path,num2str(this.currJobid));
      save([savefilename 'SC.mat' ],'SCNorm','SCMetr');                         
      save([savefilename 'SCh.mat'],'hSC','hMetr');                         % hSC is row-normalized
      
      
      this.params.ConnectomeSim.outFilenames = 'results';
            
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

