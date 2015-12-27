classdef ConnectomeAnalysis < Gridjob
  %ConnectomeAnalysis Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ConnectomeAnalysis(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ConnectomeAnalysis.pre.ConnectomeMetrics = true;
      this.params.ConnectomeAnalysis.pre.ConnectomeSim = true;
      this.params.ConnectomeAnalysis.pre.CompareWithEEG_ConnFC = true;
      this.params.ConnectomeAnalysis.pre.CompareWithEEG = true;
      
      this.params.ConnectomeAnalysis.loadSp = 0.6;                          
      this.params.ConnectomeAnalysis.loadHeur = 1;                         
      %this.params.ConnectomeAnalysis.loadHscale = 1;
      %this.params.ConnectomeAnalysis.loadKscale = 1;  
      this.params.ConnectomeAnalysis.graphH0 = 1000;
      
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
      
      param = this.params.ConnectomeAnalysis;
      
      paths = dataPaths();
      
      % comparisons over prediction performance
      if param.pre.CompareWithEEG == true
        load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics', 'CompareWithEEG', 'compareSimExp5.mat'))
        %sim.spec.dimLabels        
        
        path = strcat(paths.localTempDir,'/Results/','sp',num2str(param.loadSp),'/heur',num2str(param.loadHeur));
        if ~exist(path,'dir')
          mkdir(path);
        end
           
        % check varied parameters and access corresponding dimensions
%         dimParam = ones(4,1);
%         dimParam(1,1) = findStrcmp(overFreq.spec.dimName,'loadSp');
%         dimParam(2,1) = findStrcmp(overFreq.spec.dimName,'loadHeur');
%         dimParam(3,1) = findStrcmp(overFreq.spec.dimName,'loadHscale');
%         dimParam(4,1) = findStrcmp(overFreq.spec.dimName,'k');
% 
%         [A,B] = sort(dimParam);
%         
%         dimParam(B(A>0),:)
%         overFreq.spec.dimSize
%        
%         thisPar = [uint8(param.loadSp*10+1), param.loadHeur];
        
        if(length(overFreq.spec.dimName) == 4)                              %  Sparse x Heur x hScale x kScale
          abc=squeeze(overFreq.coh.rho(uint8(param.loadSp*10+1),param.loadHeur,:,:)); 
        end
        
        if(length(overFreq.spec.dimName) == 3)                              % Sparse x hScale x kScale, only 1 Heur
          abc=squeeze(overFreq.coh.rho(uint8(param.loadSp*10+1),:,:));      %%%% HAS TO BE FIXED
                                                                            %%%% DOES NOT WORK WHEN NOT STARTING FROM sp=0
                                                                            %%%% e.g. cannot index 7
          % strcmp(overFreq.spec.dimName(1), 'loadSp')
          % strcmp(overFreq.spec.dimName(2), 'loadHscale')
          % strcmp(overFreq.spec.dimName(3), 'k')
        end
                                                                           
        [maxcorr,idx] = max(abc(:));  
        [h,k] = ind2sub(size(abc),idx);
        imagesc(abc)                                                        % plot hScale vs kScale
        colorbar()
        title(strcat('max.corr.:',num2str(maxcorr),' for h=',num2str(h),', k=',num2str(k)))
        print(strcat(path,'/perfHvsK'),'-dpng')
        % ConnFC: ROI x ROI x Freq x Sparse x Heur x hScale x kScale      
      end
      
      % comparisons over graph metrics -- the metrics are taken BEFORE graph normalization!
      if param.pre.ConnectomeMetrics == true
        load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics', 'temp_ConnectomeMetrics', 'jobDesc.mat'))
        % load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics', 'ConnectomeMetrics', 'lotsOfFiles'))
        
        if param.graphH0 > 0
          % compare over small worldness (norm shortest paths, norm clustering coeff) etc.
          % plot normCCoeff and normShortestPath as functions of h, comparable to Strogatz & Watts 1997
        end
        
        % compare over correlation with distance matrix: hMetr.perConn.euclDist
        
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

