classdef Load3Ddataset < Gridjob
  %PREPROCESS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
%     params %possible to redefine?
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = Load3Ddataset(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% EDIT HERE: standard parameters for the job %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this.params.Load3Ddataset.outActFolder = 'inputImages'; %relative to the workpath
      this.params.Load3Ddataset.useMasks = false;
      this.params.Load3Ddataset.cutOversized = false;
      this.params.Load3Ddataset.fillUndersized = true;
      
      
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
      
      params = this.params.Load3Ddataset;
      
      data = dataPaths( );
      dbPath = fullfile(data.databases,'3Ddataset');
      dbInfo = load(fullfile(dbPath,'dbInfo.mat'));
      
      for k=1:length(dbInfo.relPaths)
        disp(k)
        currPath = dbInfo.relPaths{k};
        img = imread(fullfile(dbPath,currPath{:}));
        img = im2double(img);
        
        if params.useMasks
          maskPath = currPath(1:end-1);
          maskPath{end+1} = 'mask';
          maskPath{end+1} = [currPath{end}(1:end-4) '.mask'];
          [Data, Size] = ReadPointsData(fullfile(dbPath,maskPath{:}));
          img = img.*repmat(Data/255,[1,1,3]);
        end
        
        if params.cutOversized
          if ~isequal(size(img),[300 400 3])
            aspect = size(img,1) / size(img,2);
            if aspect > 3/4
              %% second Dim is too small -> fit second Dim
              scaleFactor = 400/size(img,2);
              img = imresize(img, [round(scaleFactor*size(img,1)) 400],'nearest');
              cutTotal = size(img,1)-300;
              cut1 = floor(cutTotal/2);
              cut2 = cutTotal - cut1;
              img = img(1+cut1:end-cut2,:,:);
            else
              %% first Dim is too small -> fit first Dim
              scaleFactor = 300/size(img,1);
              img = imresize(img, [300 round(scaleFactor*size(img,2))],'nearest');
              cutTotal = size(img,2)-400;
              cut1 = floor(cutTotal/2);
              cut2 = cutTotal - cut1;
              img = img(:,1+cut1:end-cut2,:);
            end
            
          end
          disp(size(img));
        end
        
        if params.fillUndersized
          if ~isequal(size(img),[300 400 3])
            aspect = size(img,1) / size(img,2);
            if aspect > 3/4
              %% second Dim is too small -> scale first Dim to fit and fill second Dim with gray
              scaleFactor = 300/size(img,1);
              imgRescaled = imresize(img, [300 round(scaleFactor*size(img,2))]);
              fillTotal = 400 - size(imgRescaled,2);
              img = 0.5*ones(300,400,3);
              img(:,1+floor(fillTotal/2)+(1:size(imgRescaled,2)),:) = imgRescaled;
            else
              %% first Dim is too small -> scale second Dim to fit and fill first Dim with gray
              scaleFactor = 400/size(img,2);
              imgRescaled = imresize(img, [round(scaleFactor*size(img,1)) 400]);
              fillTotal = 300 - size(imgRescaled,1);
              img = 0.5*ones(300,400,3);
              img(1+floor(fillTotal/2)+(1:size(imgRescaled,1)),:,:) = imgRescaled;
            end
            
            if min(img(:))<0
              img(img<0) = 0;
            end
            if max(img(:))>1
              img(img>1) = 1;
            end
            
            disp([min(img(:)) max(img(:))]);
            disp(size(img));
          end
        end
        
        act = img;
        outfolder = fullfile(this.workpath,this.params.Load3Ddataset.outActFolder, currPath{1:end-1});
        mkdir(outfolder);
        save(fullfile(outfolder,[currPath{end}(1:end-4) '.mat']),'act');
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

