classdef PlotRfImagePatches < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = PlotRfImagePatches(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.PlotRfImagePatches.inActImageFolder = 'labelMeInput'; %relative to the workpath
      this.params.PlotRfImagePatches.inActImageFilenames = 'act.*.mat';
      this.params.PlotRfImagePatches.inActFolder = 'layer2act'; %relative to the workpath
      this.params.PlotRfImagePatches.inActFilenames = 'act.*.mat';
      this.params.PlotRfImagePatches.inActFeatureIds = [];
      this.params.PlotRfImagePatches.plotNumPatches = 10;
      this.params.PlotRfImagePatches.plotPatchHalfLength = 12; %measured in units of the image layer
      this.params.PlotRfImagePatches.excludePatchHalfLength = 10; %measured in units of the activity layer
      this.params.PlotRfImagePatches.outPatchesFolder = 'rfImagePatches';
      this.params.PlotRfImagePatches.plotBorder = 2;

      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
            
      %list available activity files which match the given string
      inputfolder = fullfile(this.workpath,this.params.PlotRfImagePatches.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.PlotRfImagePatches.inActFilenames );
      
      inputfolderImage = fullfile(this.workpath,this.params.PlotRfImagePatches.inActImageFolder);
      [ pathlistImage, filelistImage ] = dirrec( inputfolderImage,this.params.PlotRfImagePatches.inActImageFilenames );
      
      excl = this.params.PlotRfImagePatches.excludePatchHalfLength;
      incl = this.params.PlotRfImagePatches.plotPatchHalfLength;
      
      
      
      for fileid = 1:length(filelist)
        disp(fullfile(pathlistImage{fileid},filelistImage{fileid}));
        
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        act = act.act;
        
        actImgTmp = load(fullfile(pathlistImage{fileid},filelistImage{fileid}));
        ratioX = size(actImgTmp.act,2) / size(act,2);
        ratioY = size(actImgTmp.act,1) / size(act,1);
        
        actImg = zeros(size(actImgTmp.act)+[2*incl 2*incl 0]);
        actImg(incl+1:incl+size(actImgTmp.act,1), incl+1:incl+size(actImgTmp.act,2), :) = actImgTmp.act;
        
        
        
        if fileid==1
          if isempty(this.params.PlotRfImagePatches.inActFeatureIds)
            inActFeatureIds = 1:size(act,3);
          else
            inActFeatureIds = this.params.PlotRfImagePatches.inActFeatureIds;
          end
          patches = cell(length(inActFeatureIds),1);
          for fid=1:length(inActFeatureIds)
            patches{fid} = cell(this.params.PlotRfImagePatches.plotNumPatches,1);
            for pid=1:length(patches{fid})
              patches{fid}{pid}.actValue = 0;
              patches{fid}{pid}.imgPatch = [];
              patches{fid}{pid}.imgId = 0;
              patches{fid}{pid}.imgPosX = 0;
              patches{fid}{pid}.imgPosY = 0;
            end
          end
        end
        
        for fid=1:length(inActFeatureIds)
          featureMap = act(:,:,inActFeatureIds(fid));
          
          for pid = 1:this.params.PlotRfImagePatches.plotNumPatches
            
            [maxAct,I] = max(featureMap(:));
            if maxAct > patches{fid}{pid}.actValue % better candidate found
              [y,x] = ind2sub(size(featureMap),I);
              imgx = round(x * ratioX);
              imgy = round(y * ratioY);
              
              % exclude region for continued search:
              featureMap(max(1,y-excl):min(size(featureMap,1),y+excl),...
                         max(1,x-excl):min(size(featureMap,2),x+excl)) = 0;
              
              % move all found patches with lower act:
              patches{fid}(pid+1:end) = patches{fid}(pid:end-1);
              
              % add new patch:
              patches{fid}{pid}.imgPosY = imgy;
              patches{fid}{pid}.imgPosX = imgx;
              try
                patches{fid}{pid}.imgPatch = actImg(imgy:imgy+2*incl, imgx:imgx+2*incl, :);
              catch
                pause(1);
              end
              patches{fid}{pid}.imgId = fileid;
              patches{fid}{pid}.actValue = maxAct;
            end
          end
          
        end
        
        
      end
      
      imagepatches = cellfun( @(y) cellfun(@(x) x.imgPatch, y, 'UniformOutput', false), patches, 'UniformOutput', false)';
      imagepatches = cat(2, imagepatches{:});
      imagedata = cellfun(@(x) cat(1,x,zeros(this.params.PlotRfImagePatches.plotBorder,size(x,2),size(x,3))),imagepatches, 'UniformOutput', false);
      imagedata = cellfun(@(x) cat(2,x,zeros(size(x,1),this.params.PlotRfImagePatches.plotBorder,size(x,3))),imagedata, 'UniformOutput', false);
      imagedata = cell2mat(imagedata);
      
      savedir = fullfile(this.workpath,this.params.PlotRfImagePatches.outPatchesFolder);
      if this.numJobs > 1
        savedir = fullfile(savedir,num2str(this.currJobid));
      end
      mkdir(savedir);
      imwrite(imagedata,fullfile(savedir,'patches.png'),'png')
      save(fullfile(savedir,'patches.mat'),'patches','imagepatches');
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

