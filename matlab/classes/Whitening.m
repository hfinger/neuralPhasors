classdef Whitening < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = Whitening(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.Whitening.inActFolder = 'labelMeInput'; %relative to the workpath
%       this.params.Whitening.outActFolder = 'whitenedData'; %relative to the workpath
      this.params.Whitening.inActFilenames = 'act.*.mat';
      this.params.Whitening.outWeightsFolder = 'whitening'; %relative to the workpath
      this.params.Whitening.inNumChannels = 3;
      this.params.Whitening.saveSingularValues = false;
      this.params.Whitening.maxCorrLengthDim1 = 32;
      this.params.Whitening.maxCorrLengthDim2 = 32;
      this.params.Whitening.convKernelDim1 = 21;
      this.params.Whitening.convKernelDim2 = 21;
      this.params.Whitening.reduceToConv = true;
      this.params.Whitening.epsilon = 1e-5;
      this.params.Whitening.numPatches = 10000;
      this.params.Whitening.borderBuffer = 0; %exclude these pixels around the image
      
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
      inputfolder = fullfile(this.workpath,this.params.Whitening.inActFolder);
      [ pathlist, filelist ] = dirrec( inputfolder,this.params.Whitening.inActFilenames );
      
      
      patches = loadRandomPatches( ...
        this.params.Whitening.maxCorrLengthDim1, ...
        this.params.Whitening.maxCorrLengthDim2,...
        this.params.Whitening.inNumChannels, ...
        this.params.Whitening.numPatches, ...
        pathlist, ...
        filelist, ...
        this.params.Whitening.borderBuffer );
      
      disp('now compute the mean...')
      if this.params.Whitening.reduceToConv
        meanInput = sum(sum(sum(patches,1),2),4) / (numel(patches)/size(patches,3));
      else
        meanInput = mean(patches, 4);
      end
      
      %subtract mean:
      patches = bsxfun(@minus, patches, meanInput);
      
      %reshape to [dataDims,sample] for whitening
      patches = reshape(patches,[size(patches,1)*size(patches,2)*size(patches,3) size(patches,4)]);
      
      disp('now compute the covariance matrix...')
      sigma = patches * patches' / this.params.Whitening.numPatches;
      
      disp('now compute the eigenvalues for whitening...')
      [u, s] = svd(sigma);
      
      savedir = fullfile(this.workpath,this.params.Whitening.outWeightsFolder);
      if this.numJobs > 1
        savedir = fullfile(savedir,num2str(this.currJobid));
      end
      mkdir(savedir);
      if this.params.Whitening.saveSingularValues
        save(fullfile(savedir,'us.mat'),'u','s');
      end
      
      ZCAWhite = u * diag(1 ./ sqrt(diag(s) + this.params.Whitening.epsilon)) * u';
      
      W = ZCAWhite';
      % W has now 2 dimensions: [fin,fout]
      W = reshape(W,[...
        this.params.Whitening.maxCorrLengthDim1 ...
        this.params.Whitening.maxCorrLengthDim2 ...
        this.params.Whitening.inNumChannels ...
        this.params.Whitening.maxCorrLengthDim1 ...
        this.params.Whitening.maxCorrLengthDim2 ...
        this.params.Whitening.inNumChannels]);
      % W has now 6 dimensions: [xin,yin,fin,xout,yout,fout]
      
      clear u s;
      
      
      
      % Reduce W to a convolution kernel:
      if this.params.Whitening.reduceToConv
        disp('start to reduce to a convolutional filter')
        Wconv = zeros(this.params.Whitening.convKernelDim1,this.params.Whitening.convKernelDim2,this.params.Whitening.inNumChannels,1,1,this.params.Whitening.inNumChannels);
        for i=1:this.params.Whitening.maxCorrLengthDim1-this.params.Whitening.convKernelDim1+1
          for j=1:this.params.Whitening.maxCorrLengthDim2-this.params.Whitening.convKernelDim2+1
            Wconv = Wconv + W(...
              i:i+this.params.Whitening.convKernelDim1-1,...
              j:j+this.params.Whitening.convKernelDim2-1,...
              :,...
              i+floor(this.params.Whitening.convKernelDim1/2),...
              j+floor(this.params.Whitening.convKernelDim2/2),...
              :);
          end
        end
        Wconv = Wconv / (...
          (this.params.Whitening.maxCorrLengthDim1-this.params.Whitening.convKernelDim1+1)*...
          (this.params.Whitening.maxCorrLengthDim2-this.params.Whitening.convKernelDim2+1));
        W = Wconv;
        
        conn.inputSubtract = meanInput;
        conn.shiftOutputdims = false;
        conn.W = W;
        
        save(fullfile(savedir,'weights.mat'),'-struct','conn');
      else
        
        conn.inputSubtract = meanInput;
        conn.shiftOutputdims = false;
        conn.W = W;
        
        save(fullfile(savedir,'weights.mat'),'-struct','conn');
      end
      
%       % Now whiten the inputs:
%       if ~isempty(this.params.Whitening.outActFolder)
%         disp('start to whiten the input data');
%         savedir = fullfile(this.workpath,this.params.Whitening.outActFolder);
%         if this.numJobs > 1
%           savedir = fullfile(savedir,num2str(this.currJobid));
%         end
%         mkdir(savedir);
%         
%         %list available activity files which match the given string
%         inputfolder = fullfile(this.workpath,this.params.Whitening.inActFolder);
%         [ pathlist, filelist ] = dirrec( inputfolder,this.params.Whitening.inActFilenames );
%         
%         %export whitened data to outActFolder
%         for fileid=1:length(filelist)
%           act = load(fullfile(pathlist{fileid},filelist{fileid}));
%           act = act.act;
%           if iscell(act)
%             for cellid=1:numel(act)
%               act{cellid} = applyWhitening(this, act{cellid}, W, meanInput);
%             end
%           else
%             act = applyWhitening(this, act, W, meanInput);
%           end
%           savepath = fullfile(savedir, pathlist{fileid}(length(inputfolder)+2:end));
%           mkdir(savepath);
%           save(fullfile(savepath,filelist{fileid}),'act');
%         end
%       end
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
%     function act = applyWhitening(this, act, W, meanInput)
%       if this.params.Whitening.reduceToConv
%         act = bsxfun(@minus,act,meanInput);
%         act = conv3d(permute(W,[1 2 3 6 4 5]),act);
%       else
%         act = bsxfun(@minus,act,meanInput);
%         act = reshape( reshape(W,[size(W,1)*size(W,2)*size(W,3) size(W,4)*size(W,5)*size(W,6)])'*act(:) ,size(act));
%       end
%     end
    
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

