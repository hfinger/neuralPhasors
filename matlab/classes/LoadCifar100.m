classdef LoadCifar100 < Gridjob
  %PREPROCESS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
%     params %possible to redefine?
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = LoadCifar100(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% EDIT HERE: standard parameters for the job %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this.params.LoadCifar100.useTestset = true;
      this.params.LoadCifar100.imageids = []; %or [] for all
      this.params.LoadCifar100.outActFolder = 'cifarInput'; %relative to the workpath
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% EDIT HERE: laod image and resize and save %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      datapaths = dataPaths( );
      if this.params.LoadCifar100.useTestset
        data = load(fullfile(datapaths.cifar100matlab, 'test.mat'));
      else
        data = load(fullfile(datapaths.cifar100matlab, 'train.mat'));
      end
      data = data.data';
      if ~isempty(this.params.LoadCifar100.imageids)
        data = data(:,this.params.LoadCifar100.imageids);
      end
      data = reshape(data,[32 32 3 size(data,2)]);
      data = im2double(data);
      act = squeeze(num2cell(data,[1 2 3])); %#ok<NASGU>
      outfolder = fullfile(this.workpath,this.params.LoadCifar100.outActFolder);
      mkdir(outfolder);
      save(fullfile(outfolder,'act1.mat'),'act');
      
      disp('finished loading cifarImages')
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% EDIT HERE: do some clean up and saving %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

