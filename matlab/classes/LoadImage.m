classdef LoadImage < Gridjob
  %PREPROCESS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
%     params %possible to redefine?
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = LoadImage(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% EDIT HERE: standard parameters for the job %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      this.params.LoadImage.resizeToX = [];
      this.params.LoadImage.resizeToY = [];
      this.params.LoadImage.doResize = false;
      this.params.LoadImage.filepath = 'img.png';
      this.params.LoadImage.outActFolder = 'imageInput'; %relative to the workpath
      this.params.LoadImage.useFilenameAsDir = false;
      
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
      
        imgPath = this.params.LoadImage.filepath;
        act = imread(imgPath);
        act = im2double(act);
        if this.params.LoadImage.doResize
          act = imresize(act, [this.params.LoadImage.resizeToY this.params.LoadImage.resizeToX]); %#ok<NASGU>
        end
        
        if this.params.LoadImage.useFilenameAsDir
          [~, name, ext] = fileparts(imgPath);
          outfolder = fullfile(this.workpath,this.params.LoadImage.outActFolder, [name ext]);
        else
          if this.numJobs==1
            outfolder = fullfile(this.workpath,this.params.LoadImage.outActFolder);
          else
            outfolder = fullfile(this.workpath,this.params.LoadImage.outActFolder, num2str(this.currJobid));
          end
        end
        
        
        mkdir(outfolder);
        save(fullfile(outfolder,'act1.mat'),'act');
      
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

