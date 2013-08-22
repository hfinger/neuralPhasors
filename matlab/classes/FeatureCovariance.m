classdef FeatureCovariance < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = FeatureCovariance(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.FeatureCovariance.inActFolder = 'labelMeInput'; %relative to the workpath
      this.params.FeatureCovariance.inActFilenames = 'act.*.mat';
      this.params.FeatureCovariance.fileid = [];
      this.params.FeatureCovariance.outCovFolder = 'FeatureCorrCov'; %relative to the workpath
      this.params.FeatureCovariance.maxCovLengthDim1 = 32;
      this.params.FeatureCovariance.maxCovLengthDim2 = 32;
      this.params.FeatureCovariance.numSamplesPerImage = 1000;
      this.params.FeatureCovariance.borderBuffer = 0; %exclude these pixels around the image
      this.params.FeatureCovariance.saveCorr = true;
      this.params.FeatureCovariance.saveCov = true;
      this.params.FeatureCovariance.saveVar = true;
      this.params.FeatureCovariance.saveCorrPvalue = true;
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      covParam = this.params.FeatureCovariance;
      
      savepath = fullfile(this.workpath,covParam.outCovFolder);
      if this.numJobs > 1
        savepath = fullfile(savepath,num2str(this.currJobid));
      end
      mkdir(savepath);
      
      inputfolder = fullfile(this.workpath,covParam.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,covParam.inActFilenames );
      
      if ~isempty(covParam.fileid)
        filelist = filelist(covParam.fileid);
        pathlist = pathlist(covParam.fileid);
      end
      
      %% check number of features:
      act = load(fullfile(pathlist{1},filelist{1}));
      if iscell(act.act)
        inNumFeatures = size(act.act{1},3);
      else
        inNumFeatures = size(act.act,3);
      end
      
      %% First start to select random patches and compute the mean of all patches:
      patchCenters = cell(size(filelist));
      patchMean = zeros(2*covParam.maxCovLengthDim1+1,2*covParam.maxCovLengthDim2+1,inNumFeatures);
      patchMeanCenter = zeros(1,1,inNumFeatures);
      numPatches = 0;
      
      borderX = covParam.borderBuffer + covParam.maxCovLengthDim1;
      borderY = covParam.borderBuffer + covParam.maxCovLengthDim2;
      
      disp('start calculating means')
      for fileid=1:numel(filelist)
        disp(['fileid=' num2str(fileid)])
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        
        if iscell(act.act)
          actCell = act.act;
        else
          actCell{1} = act.act;
        end
        
        patchCenters{fileid} = cell(size(actCell));
        for cellid=1:length(actCell)
          patchCenters{fileid}{cellid}.x = borderX + randi(size(actCell{cellid},1)-2*borderX,covParam.numSamplesPerImage,1);
          patchCenters{fileid}{cellid}.y = borderY + randi(size(actCell{cellid},2)-2*borderY,covParam.numSamplesPerImage,1);
          
          for i=1:covParam.numSamplesPerImage
            patchMean = patchMean + actCell{cellid}(...
              patchCenters{fileid}{cellid}.x(i)-covParam.maxCovLengthDim1:patchCenters{fileid}{cellid}.x(i)+covParam.maxCovLengthDim1,...
              patchCenters{fileid}{cellid}.y(i)-covParam.maxCovLengthDim2:patchCenters{fileid}{cellid}.y(i)+covParam.maxCovLengthDim2,...
              :);
            patchMeanCenter = patchMeanCenter + actCell{cellid}(...
              patchCenters{fileid}{cellid}.x(i),...
              patchCenters{fileid}{cellid}.y(i),...
              :);
          end
          
          numPatches = numPatches + covParam.numSamplesPerImage;
        end
        
      end
      
      patchMean = patchMean / numPatches;
      patchMeanCenter = patchMeanCenter / numPatches;
      save(fullfile(savepath,'patchMeans.mat'),'patchCenters','patchMean','patchMeanCenter');
      
      %% Now start to compute the covariance matrix:
      
      patchCorrCov = zeros((2*covParam.maxCovLengthDim1+1)*(2*covParam.maxCovLengthDim2+1)*inNumFeatures,inNumFeatures);
      varA = zeros(2*covParam.maxCovLengthDim1+1,2*covParam.maxCovLengthDim2+1,inNumFeatures);
      varB = zeros(1,1,inNumFeatures);
            
      disp('start calculating cov')
      for fileid=1:numel(filelist)
        disp(['fileid=' num2str(fileid)])
        act = load(fullfile(pathlist{fileid},filelist{fileid}));
        
        if iscell(act)
          actCell = act.act;
        else
          actCell{1} = act.act;
        end
        
        for cellid=1:length(actCell)
          for i=1:covParam.numSamplesPerImage
            A = actCell{cellid}(...
              patchCenters{fileid}{cellid}.x(i)-covParam.maxCovLengthDim1:patchCenters{fileid}{cellid}.x(i)+covParam.maxCovLengthDim1,...
              patchCenters{fileid}{cellid}.y(i)-covParam.maxCovLengthDim2:patchCenters{fileid}{cellid}.y(i)+covParam.maxCovLengthDim2,...
              :) - patchMean;
            B = actCell{cellid}(...
              patchCenters{fileid}{cellid}.x(i),...
              patchCenters{fileid}{cellid}.y(i),...
              :) - patchMeanCenter;
            patchCorrCov = patchCorrCov + A(:) * B(:)';
            varA = varA + A.^2;
            varB = varB + B.^2;
          end
        end
      end
      
      patchCorrCov = patchCorrCov / numPatches;
      varA = varA / numPatches;
      varB = varB / numPatches;
      stdA = sqrt(varA);
      stdB = sqrt(varB);
      
      if this.params.FeatureCovariance.saveVar
        save(fullfile(savepath,'patchVar.mat'),'-v7.3','varA','varB','stdA','stdB');
      end
      
      patchCorrCov = reshape(patchCorrCov,[2*covParam.maxCovLengthDim1+1,...
        2*covParam.maxCovLengthDim2+1,...
        inNumFeatures,...
        inNumFeatures]);
      
      if this.params.FeatureCovariance.saveCov
        save(fullfile(savepath,'patchCov.mat'),'-v7.3','patchCorrCov');
      end
      
      if this.params.FeatureCovariance.saveCorr
        disp('now compute the correlation...')
        patchCorrCov = patchCorrCov ./ bsxfun(@times,stdA,shiftdim(stdB,-1));
        
        if this.params.FeatureCovariance.saveCorrPvalue
          p = this.pvalPearson('b', patchCorrCov, numPatches);
          save(fullfile(savepath,'patchCorr.mat'),'-v7.3','patchCorrCov','p');
        else
          save(fullfile(savepath,'patchCorr.mat'),'-v7.3','patchCorrCov');
        end
      end
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function p = pvalPearson(this, tail, rho, n)
      %PVALPEARSON Tail probability for Pearson's linear correlation.
      t = rho.*sqrt((n-2)./(1-rho.^2));
      switch tail
        case 'b' % 'both or 'ne'
          p = 2*tcdf(-abs(t),n-2);
        case 'r' % 'right' or 'gt'
          p = tcdf(-t,n-2);
        case 'l' % 'left' or 'lt'
          p = tcdf(t,n-2);
      end
    end
    
    function plotCov(this)
      
      patchCorrCov(:,18,:,:)=NaN;
      patchCorrCov(18,:,:,:)=NaN;
      imagesc(reshape(permute(patchCorrCov,[1 3 2 4]),[18*108 18*108])); set(gca,'clim',[-0.2 0.5]); axis equal
      
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

