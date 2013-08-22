classdef EvaluateFixations < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = EvaluateFixations(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.EvaluateFixations.inPhaseFeatFolder = 'inPhaseFeat'; %relative to the workpath
      this.params.EvaluateFixations.inPhaseFeatFilenames = 'act.*.mat';
      this.params.EvaluateFixations.fileid = [];
      this.params.EvaluateFixations.inFixationMat = [];
      this.params.EvaluateFixations.inFixationScaleX = 1;
      this.params.EvaluateFixations.inFixationScaleY = 1;
      this.params.EvaluateFixations.inFixationOffsetX = 0;
      this.params.EvaluateFixations.inFixationOffsetY = 0;
      this.params.EvaluateFixations.inFixationFileid = [];
      this.params.EvaluateFixations.outFolder = 'outFixationEval'; %relative to the workpath
      this.params.EvaluateFixations.figureOutFolder = []; %relative to the workpath
      this.params.EvaluateFixations.compareAll2All = false;
      
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
      inputfolder = fullfile(this.workpath,this.params.EvaluateFixations.inPhaseFeatFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.EvaluateFixations.inPhaseFeatFilenames );
      
      if ~isempty(this.params.EvaluateFixations.fileid)
        filelist = filelist(this.params.EvaluateFixations.fileid);
        pathlist = pathlist(this.params.EvaluateFixations.fileid);
      end
      
      if ~isempty(this.params.EvaluateFixations.inFixationFileid)
        imageIds = this.params.EvaluateFixations.inFixationFileid;
      else
        imageIds = unique(this.params.EvaluateFixations.inFixationMat.image);
        imageIds = sort(imageIds);
      end
      
      if this.params.EvaluateFixations.compareAll2All 
        meanFixationPhaseFeat = zeros(length(filelist),length(imageIds));
      else
        meanFixationPhaseFeat = zeros(1,length(filelist));
      end
      
      tmpSum = 0;
      tmpNum = 0;
      for fileid=1:length(filelist)
        disp(['fixImageId:' num2str(imageIds(fileid)) ' ' fullfile(pathlist{fileid},filelist{fileid})]);
        
        phaseFeat = load(fullfile(pathlist{fileid},filelist{fileid}));
        phaseFeat = phaseFeat.feat;
          
        for fixationfileid=1:length(imageIds)
        
          useIds = find(this.params.EvaluateFixations.inFixationMat.image == imageIds(fixationfileid) & this.params.EvaluateFixations.inFixationMat.fix>2);

          x = (double(this.params.EvaluateFixations.inFixationMat.x(useIds)) - this.params.EvaluateFixations.inFixationOffsetX) * this.params.EvaluateFixations.inFixationScaleX;
          y = (double(this.params.EvaluateFixations.inFixationMat.y(useIds)) - this.params.EvaluateFixations.inFixationOffsetY) * this.params.EvaluateFixations.inFixationScaleY;
          
          x = round(x);
          y = round(y);

          outliers = x<=0 | y<=0 | x>size(phaseFeat,2) | y>size(phaseFeat,1);
          x(outliers) = [];
          y(outliers) = [];

          linInd = sub2ind(size(phaseFeat),y,x);
          disp(num2str(length(linInd)));
          meanFixationPhaseFeat(fileid,fixationfileid) = mean(phaseFeat(linInd));
          
          if fixationfileid==fileid
            tmpSum = tmpSum + sum(phaseFeat(linInd));
            tmpNum = tmpNum + length(phaseFeat(linInd));
          end

          if ~isempty(this.params.EvaluateFixations.figureOutFolder)
            if this.numJobs > 1
              savepathFigure = fullfile(this.workpath,this.params.EvaluateFixations.figureOutFolder,num2str(this.currJobid));
            else
              savepathFigure = fullfile(this.workpath,this.params.EvaluateFixations.figureOutFolder);
            end
            imagesc(phaseFeat);
            hold on;
            plot(x,y,'LineWidth',4,'Color','r','Marker','o','LineStyle','none')
            hold off;
            mkdir(savepathFigure);
            print(gcf,'-dpsc2',fullfile(savepathFigure,['fixations' num2str(fileid) '.eps']));
          end
        end
        
      end
      
      tmpSum = tmpSum / tmpNum;
      
      disp(num2str(tmpNum));
      disp(num2str(tmpSum));
      
      if this.numJobs > 1
        savepath = fullfile(this.workpath,this.params.EvaluateFixations.outFolder,num2str(this.currJobid));
      else
        savepath = fullfile(this.workpath,this.params.EvaluateFixations.outFolder);
      end
      
      mkdir(savepath);
      save(fullfile(savepath),'-v7.3','meanFixationPhaseFeat');
      
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

