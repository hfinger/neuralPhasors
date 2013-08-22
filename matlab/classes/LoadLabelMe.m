classdef LoadLabelMe < Gridjob
  %PREPROCESS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
%     params %possible to redefine?
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = LoadLabelMe(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% EDIT HERE: standard parameters for the job %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      this.params.LoadLabelMe.resizeToX = 400;
      this.params.LoadLabelMe.resizeToY = 300;
      this.params.LoadLabelMe.doResize = true;
      this.params.LoadLabelMe.interpolateNearestNeighbor = false; %set to true to constrain the output range between 0 and 1
      this.params.LoadLabelMe.catName = '05june05_static_street_porter';
      this.params.LoadLabelMe.fileid = 1; %or [] for all in this category or a vector with ids
      this.params.LoadLabelMe.outActFolder = 'labelMeInput'; %relative to the workpath
      this.params.LoadLabelMe.convertToGray = false;
      
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
      
      data = dataPaths( );
      files = dir(fullfile(data.HOMEIMAGES, this.params.LoadLabelMe.catName, '*.jpg'));
      files = {files.name};
      if ~isempty(this.params.LoadLabelMe.fileid)
        files=files(this.params.LoadLabelMe.fileid);
      end
      for j=1:length(files)
        imgPath = fullfile(data.HOMEIMAGES, this.params.LoadLabelMe.catName, files{j});
        act = imread(imgPath);
        if 0.75~=size(act,1)/size(act,2)
          disp(['WARNING: Apect Ratio of Image is Changed: ' imgPath]);
        end
        if this.params.LoadLabelMe.convertToGray
          act = rgb2gray(act);
        end
        act = im2double(act);
        if this.params.LoadLabelMe.doResize
          if this.params.LoadLabelMe.interpolateNearestNeighbor
            act = imresize(act, [this.params.LoadLabelMe.resizeToY this.params.LoadLabelMe.resizeToX],'nearest'); %#ok<NASGU>
          else
            act = imresize(act, [this.params.LoadLabelMe.resizeToY this.params.LoadLabelMe.resizeToX]); %#ok<NASGU>
          end
        end
        outfolder = fullfile(this.workpath,this.params.LoadLabelMe.outActFolder, this.params.LoadLabelMe.catName, files{j});
        mkdir(outfolder);
        save(fullfile(outfolder,'act1.mat'),'act');
      end
      disp('finished loading labelMeImages')
      
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

