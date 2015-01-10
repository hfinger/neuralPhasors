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
      this.params.ConnectomeAnalysis.loadHscale = 1;
      this.params.ConnectomeAnalysis.loadKscale = 1;      
      
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
      
      if param.pre.CompareWithEEG == true
        load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics', 'CompareWithEEG', 'compareSimExp5.mat'))
        %sim.spec.dimLabels
        
        path = strcat(paths.localTempDir,'/Results/','sp',num2str(param.loadSp),'/heur',num2str(param.loadHeur));
        if ~exist(path,'dir')
          mkdir(path);
        end
        
        abc=squeeze(overFreq.coh.rho(uint8(param.loadSp*10+1),param.loadHeur,:,:)); % Sparse x Heur x hScale x kScale
        imagesc(abc)                                                         % plot hScale vs kScale
        colorbar()
        print(strcat(path,'/perfHvsK'),'-dpng')
        
        % ConnFC: ROI x ROI x Freq x Sparse x Heur x hScale x kScale
        
      end
      
%       %%%
%       path = fullfile(this.workpath,param.outFilenames);
%       if ~exist(path,'dir')
%         mkdir(path);
%       end
%       savefilename = fullfile(path,num2str(this.currJobid));
%       %save([savefilename 'SC.mat' ],'SCNorm','SCMetr');                    % removed: take case hScale = 0 as SC matrix, only when add == true in homotopHeur!
%       save([savefilename 'SC.mat'],'hSC','hMetr');                          % hSC is row-normalized
%       
%       
%       this.params.ConnectomeSim.outFilenames = 'results';
      
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

