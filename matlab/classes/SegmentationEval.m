classdef SegmentationEval < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = SegmentationEval(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.SegmentationEval.inActFolder = 'layer1Act';
      this.params.SegmentationEval.inActFilenames = 'act1.mat';
      this.params.SegmentationEval.inPhaseFolder = 'labelMePhase';
      this.params.SegmentationEval.outSegEvalFolder = 'layer1SegEval';
      this.params.SegmentationEval.catName = '05june05_static_street_boston';
      this.params.SegmentationEval.fileid = 1:50;
      
      this.params.SegmentationEval.minPixelPerSeg = 50;
      this.params.SegmentationEval.maxPixelPerSeg = Inf;
      this.params.SegmentationEval.borderSize=0;
      this.params.SegmentationEval.numPairs=100000; % set numPairs=0 to calc only probSegSync
      this.params.SegmentationEval.numSubpop=1000;
      this.params.SegmentationEval.numSamplesInSubpop=1000;
      this.params.SegmentationEval.phasePerPixel=false;
      this.params.SegmentationEval.time=20;
      this.params.SegmentationEval.initRandstream = true;
      this.params.SegmentationEval.pairOnlyWithOneNonmatching = false;
      this.params.SegmentationEval.verboseLevel = 0;
      this.params.SegmentationEval.simOtherVarParam = [];
      this.params.SegmentationEval.skipMissingImg = false;
      this.params.SegmentationEval.borderEvalNumPerImg = 10;
      this.params.SegmentationEval.borderEvalNghSize = 10;
      this.params.SegmentationEval.borderAngleSigma = 3;
      this.params.SegmentationEval.borderAngleUseActivityWeighting = true;
      this.params.SegmentationEval.borderSyncUseActivityWeighting = true;
      this.params.SegmentationEval.useHueInsteadPhase = false;
      this.params.SegmentationEval.segSizeBinEdges = exp(4:11);
      
      this.params.SegmentationEval.resizeToX = 400;
      this.params.SegmentationEval.resizeToY = 300;
      this.params.SegmentationEval.cropX = [];
      this.params.SegmentationEval.cropY = [];
      
      this.params.SegmentationEval.doCalcSegSync = true;
      this.params.SegmentationEval.doCalcBorderSync = true;
      this.params.SegmentationEval.doCalcBorderAngle = true;
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    %% Start: is executed before all individual parameter jobs are started
    function startJob(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
      
      params = this.params.SegmentationEval;
      
      data = dataPaths( );
      labelMeFolder = params.catName;
      
      files = dir(fullfile(data.HOMEIMAGES,labelMeFolder, '*.jpg'));
      files = {files.name};
      if ~isempty(params.fileid)
        files = files(params.fileid);
      end
      
      addpath(data.labelMeToolbox);
      
      %% read segmentation masks of all simulated images:
      
      % Variables for seg eval:
      xpoly = cell(length(files),1);
      ypoly = cell(length(files),1);
      segmentSize = cell(length(files),1);
      linelength = cell(length(files),1);
      
      %% start loop over masks:
      for maskfileIdx=1:length(files)
        D = LMdatabase(data.HOMEANNOTATIONS,data.HOMEIMAGES, {labelMeFolder}, files(maskfileIdx));
        tmpImg = imread(fullfile(data.HOMEIMAGES, labelMeFolder, files{maskfileIdx} ));
        imagesizeX = size(tmpImg,2);
        imagesizeY = size(tmpImg,1);
        
        xpoly{maskfileIdx} = cell(length(D.annotation.object),1);
        ypoly{maskfileIdx} = cell(length(D.annotation.object),1);
        linelength{maskfileIdx} = cell(length(D.annotation.object),1);
        segmentSize{maskfileIdx} = zeros(length(D.annotation.object),1);
        polyidx = 1;
        for obj=1:length(D.annotation.object)
          if str2double(D.annotation.object(obj).deleted)
            continue;
          end
          xtmp=double(D.annotation.object(obj).polygon.x) * params.resizeToX / imagesizeX;
          ytmp=double(D.annotation.object(obj).polygon.y) * params.resizeToY / imagesizeY;
          if ~isempty(params.cropX)
            xtmp=xtmp-params.cropX;
          end
          if ~isempty(params.cropY)
            ytmp=ytmp-params.cropY;
          end
          BW = poly2mask(xtmp, ytmp, params.resizeToY, params.resizeToX);
          SE = strel('diamond',params.borderSize);
          BW = imerode(BW,SE);
          if sum(BW(:))<params.minPixelPerSeg || sum(BW(:))>params.maxPixelPerSeg
            continue;
          end
          segmentSize{maskfileIdx}(polyidx)=sum(BW(:));
          xpoly{maskfileIdx}{polyidx} = xtmp;
          ypoly{maskfileIdx}{polyidx} = ytmp;
          
          %calc line lengths:
          linelength{maskfileIdx}{polyidx} = sqrt(diff(xpoly{maskfileIdx}{polyidx}([1:end 1])).^2+diff(ypoly{maskfileIdx}{polyidx}([1:end 1])).^2);
          
          polyidx = polyidx + 1;
        end
        xpoly{maskfileIdx}(polyidx:end) = [];
        ypoly{maskfileIdx}(polyidx:end) = [];
        linelength{maskfileIdx}(polyidx:end) = [];
        segmentSize{maskfileIdx}(polyidx:end) = [];
      end
      save(fullfile(this.temppath,'segments.mat'),'xpoly','ypoly','segmentSize','linelength','files');
      
      %%%% END EDIT HERE:                                        %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Now start the job via superclass
      %       start@Gridjob(this)
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      params = this.params.SegmentationEval;
      
      if params.initRandstream
        if isinteger(params.initRandstream)
          reset(RandStream.getGlobalStream,params.initRandstream);
        elseif ischar(params.initRandstream) && strcmp(params.initRandstream,'clock')
          reset(RandStream.getGlobalStream,sum(100*clock));
        else %is boolean
          reset(RandStream.getGlobalStream,this.currJobid);
        end
      end
      
      tmp = load(fullfile(this.temppath,'segments.mat'));
      xpoly = tmp.xpoly;
      ypoly = tmp.ypoly;
      segmentSize = tmp.segmentSize;
      linelength = tmp.linelength;
      files = tmp.files;
      
      %% select positions for border eval:
      if params.doCalcBorderSync
        % Variables for border eval:
        borders.borderPosX = cell(length(files),1);
        borders.borderPosY = cell(length(files),1);
        borders.borderPolyId = cell(length(files),1);
        borders.borderPolyLineId = cell(length(files),1);
        borders.borderPolyLineOffset = cell(length(files),1);
        borders.borderAngle = cell(length(files),1);
        for maskfileIdx=1:length(files)
          
          disp(['maskfileIdx=' num2str(maskfileIdx)])
          
          if iscell(params.borderEvalNghSize)
            maxNghSize = max(cell2mat(params.borderEvalNghSize));
          else
            maxNghSize = params.borderEvalNghSize;
          end
          
          numLinesPerPoly = cellfun(@length,linelength{maskfileIdx});
          polyLineIdEdges = [1;cumsum(numLinesPerPoly)+1];
          allLineLengths = cell2mat(linelength{maskfileIdx});
          allLineEdges = [0;cumsum(allLineLengths)];
          totalLength = allLineEdges(end);
          allLineEdges(end)=allLineEdges(end)+1;
          
          borderPosX=zeros(1,params.borderEvalNumPerImg);
          borderPosY=zeros(1,params.borderEvalNumPerImg);
          borderAngle=zeros(1,params.borderEvalNumPerImg);
          borderPolyLineId=zeros(1,params.borderEvalNumPerImg);
          borderPolyLineOffset=zeros(1,params.borderEvalNumPerImg);
          borderPolyId=zeros(1,params.borderEvalNumPerImg);
          
          id=1;
          while id<=params.borderEvalNumPerImg
            r = totalLength*rand(1,1);
            [~,borderPolyLineId(id)] = histc(r,allLineEdges);
            borderPolyLineOffset(id) = r - allLineEdges(borderPolyLineId(id));
            [~,borderPolyId(id)] = histc(borderPolyLineId(id),polyLineIdEdges);
            borderPolyLineId(id) = borderPolyLineId(id) - polyLineIdEdges(borderPolyId(id)) + 1;
            
            startVertexId=borderPolyLineId(id);
            endVertexId=startVertexId+1;
            if endVertexId>numLinesPerPoly(borderPolyId(id))
              endVertexId=1;
            end
            
            %Line endpoints:
            xcoords=xpoly{maskfileIdx}{borderPolyId(id)}([startVertexId endVertexId]);
            ycoords=ypoly{maskfileIdx}{borderPolyId(id)}([startVertexId endVertexId]);
            
            lineFraction = borderPolyLineOffset(id) / linelength{maskfileIdx}{borderPolyId(id)}(borderPolyLineId(id));
            borderPosX(id) = (1-lineFraction)*xcoords(1)+lineFraction*xcoords(2);
            borderPosY(id) = (1-lineFraction)*ycoords(1)+lineFraction*ycoords(2);
            
            borderAngle(id) = atan2(ycoords(2)-ycoords(1),xcoords(2)-xcoords(1));
            
            if borderPosX(id)>maxNghSize && ...
                borderPosY(id)>maxNghSize && ...
                borderPosX(id)+maxNghSize<params.resizeToX && ...
                borderPosY(id)+maxNghSize<params.resizeToY
              id=id+1;
            end
          end
          
          borders.borderPosX{maskfileIdx}=borderPosX;
          borders.borderPosY{maskfileIdx}=borderPosY;
          borders.borderPolyId{maskfileIdx}=borderPolyId;
          borders.borderPolyLineId{maskfileIdx}=borderPolyLineId;
          borders.borderPolyLineOffset{maskfileIdx}=borderPolyLineOffset;
          borders.borderAngle{maskfileIdx}=borderAngle;
        end
      end
      
      %%
      numsims = length(this.params.SegmentationEval.fileid);
      
      if params.pairOnlyWithOneNonmatching
        disp('select paired segments and border regions...')
        pairMaskFileWithId=segmentSize;
        for maskid=1:numsims
          for objid=1:length(pairMaskFileWithId{maskid})
            pairMaskFileWithId{maskid}(objid)=randi(numsims-1); %random nonmatching
            if pairMaskFileWithId{maskid}(objid) >= maskid
              pairMaskFileWithId{maskid}(objid) = pairMaskFileWithId{maskid}(objid)+1;
            end
          end
        end
        
        if params.borderEvalNumPerImg>0
          pairBorderObjWithImgId=cell(numsims,1);
          for maskid=1:numsims
            for objid=1:params.borderEvalNumPerImg
              pairBorderObjWithImgId{maskid}(objid)=randi(numsims-1); %random nonmatching
              if pairBorderObjWithImgId{maskid}(objid) >= maskid
                pairBorderObjWithImgId{maskid}(objid) = pairBorderObjWithImgId{maskid}(objid)+1;
              end
            end
          end
        end
      end
      
      %% calc all segmentation masks and apply them to all simulations results:
      tmpMatching = 0;
      tmpNonMatching = 0;
      tmpMatchingCounter = 0;
      tmpNonMatchingCounter = 0;
      
      tmpAngleErrorMatching = 0;
      tmpAngleErrorNonMatching = 0;
      tmpAngleErrorMatchingCounter = 0;
      tmpAngleErrorNonMatchingCounter = 0;
      
      output = cell(numsims,numsims);
      outputBorders = cell(numsims,numsims);
      for simIdx=1:numsims
        disp(['simidx: ' num2str(simIdx) ' / ' num2str(numsims)]);
        
        if params.useHueInsteadPhase || params.verboseLevel>2
          % Load rgb image
          data = dataPaths( );
          tmpImg = imread(fullfile(data.HOMEIMAGES,params.catName, files{simIdx} ));
          tmpImg = imresize(tmpImg,[params.resizeToY params.resizeToX]);
          if ~isempty(params.cropX)
            tmpImg=tmpImg(:,params.cropX+1:end-params.cropX,:);
          end
          if ~isempty(params.cropY)
            tmpImg=tmpImg(params.cropY+1:end-params.cropY,:,:);
          end
        end
        
        if params.useHueInsteadPhase
          hsv_image = rgb2hsv(tmpImg);
          phase = 2*pi*hsv_image(:,:,1);
          features = ones(size(phase));
        else
          features = load(fullfile(this.workpath,params.inActFolder,params.catName, files{simIdx},'act1.mat'));
          features = features.act;
          phase = load(fullfile(this.workpath,params.inPhaseFolder,params.catName, files{simIdx},['phaseIter' num2str(params.time) '.mat']));
          phase = phase.phase;
        end
        
        if params.verboseLevel>2
          circSegVar = circ_mean(phase, features, 3);
        end
        
        if params.phasePerPixel
          phase = circ_mean(phase, features, 3);
          features = ones(size(phase));
        end
        
        for maskfileIdx=1:numsims
          tic;
          if params.doCalcSegSync
            output{simIdx,maskfileIdx} = cell(1,length(xpoly{maskfileIdx}));
            for obj=1:length(xpoly{maskfileIdx})
              
              if params.pairOnlyWithOneNonmatching && maskfileIdx~=simIdx && pairMaskFileWithId{maskfileIdx}(obj)~=simIdx
                continue;
              end
              
              if params.verboseLevel>0
                disp(['simIdx=' num2str(simIdx) ' maskfileIdx=' num2str(maskfileIdx) ' obj=' num2str(obj) ])
              end
              
              BW = poly2mask(xpoly{maskfileIdx}{obj}, ypoly{maskfileIdx}{obj}, size(features,1), size(features,2));
              %imagesc(BW)
              
              %% Now select the segment-core and the neighborhood regions:
              SEborder = strel('diamond',params.borderSize);
              BWplusBorder = imdilate(BW,SEborder);
              
              if params.borderSize>0
                BW = imerode(BW,SEborder);
              end
              
              for nghSize = 1:40
                SE = strel('diamond',nghSize);
                NH = imdilate(BWplusBorder,SE);
                NH = logical(NH - BWplusBorder);
                numPixelInNgh = sum(NH(:));
                if numPixelInNgh > sum(BW(:))
                  break;
                end
              end
              
              %% Evaluate the statistics
              output{simIdx,maskfileIdx}{obj} = this.calcSyncValues( features, phase, BW, NH );
              
              if params.verboseLevel>2
                disp(output{simIdx,maskfileIdx}{obj})
                
                figure(1);
                imagesc(tmpImg);
                
                figure(2);
                BW2=BW(:,:,1);
                circSegVar2 = zeros(size(circSegVar));
                circSegVar2(BW2(:)) = circSegVar(BW2(:));
                imagesc(circSegVar2);
                colormap hsv;
                caxis([-pi pi]);
                colorbar;
                
                figure(3);
                imagesc(circSegVar);
                colormap hsv;
                caxis([-pi pi]);
                colorbar;
                
                figure(4);
                circSegVar3 = zeros(size(circSegVar));
                circSegVar3(NH(:)) = circSegVar(NH(:));
                imagesc(circSegVar3);
                colormap hsv;
                caxis([-pi pi]);
                colorbar;
                
                drawnow;
                pause(0.1);
              end
            end
          end
          
          %% now evaluate borders:
          outputBorders{simIdx,maskfileIdx} = cell(1,params.borderEvalNumPerImg);
          
          for obj=1:params.borderEvalNumPerImg
            if params.doCalcBorderSync
              if params.pairOnlyWithOneNonmatching && maskfileIdx~=simIdx && pairBorderObjWithImgId{maskfileIdx}(obj)~=simIdx
                continue;
              end
              
              cx=borders.borderPosX{maskfileIdx}(obj);
              cy=borders.borderPosY{maskfileIdx}(obj);
              ix=size(features,2);
              iy=size(features,1);
              r=params.borderEvalNghSize;
              [x,y]=meshgrid(-(cx-1):(ix-cx),-(cy-1):(iy-cy));
              mask=((x.^2+y.^2)<=r^2 & mod(atan2(y,x)-borders.borderAngle{maskfileIdx}(obj),2*pi)<pi );
              maskFull=repmat(mask,[1 1 size(phase,3)]);
              if params.borderSyncUseActivityWeighting
                meanPhase1 = circ_mean(phase(maskFull), features(maskFull), 1);
              else
                meanPhase1 = circ_mean(phase(maskFull), [], 1);
              end
              mask=((x.^2+y.^2)<=r^2 & mod(atan2(y,x)-borders.borderAngle{maskfileIdx}(obj),2*pi)>pi );
              maskFull=repmat(mask,[1 1 size(phase,3)]);
              if params.borderSyncUseActivityWeighting
                meanPhase2 = circ_mean(phase(maskFull), features(maskFull), 1);
              else
                meanPhase2 = circ_mean(phase(maskFull), [], 1);
              end
              phaseDiff=mod(meanPhase1-meanPhase2,2*pi);
              if phaseDiff>pi
                phaseDiff=abs(phaseDiff-2*pi);
              end
              
              outputBorders{simIdx,maskfileIdx}{obj}.phaseDiff = phaseDiff;
              if maskfileIdx==simIdx
                tmpMatching = tmpMatching+phaseDiff;
                tmpMatchingCounter = tmpMatchingCounter + 1;
              else
                tmpNonMatching = tmpNonMatching+phaseDiff;
                tmpNonMatchingCounter = tmpNonMatchingCounter + 1;
              end
              
              if params.verboseLevel>2
                disp(num2str(phaseDiff));
                
                polyId = borders.borderPolyId{maskfileIdx}(obj);
                BW = poly2mask(xpoly{maskfileIdx}{polyId}, ypoly{maskfileIdx}{polyId}, size(features,1), size(features,2));
                
                figure(1);
                imagesc(BW); colormap gray;
                
                figure(2);
                imagesc(circSegVar);
                colormap hsv;
                caxis([-pi pi]);
                colorbar;
                hold on
                th = 0:pi/50:2*pi;
                xunit = params.borderEvalNghSize * cos(th) + borders.borderPosX{maskfileIdx}(obj);
                yunit = params.borderEvalNghSize * sin(th) + borders.borderPosY{maskfileIdx}(obj);
                plot(xunit, yunit,'o');
                hold off
                
                figure(3)
                imagesc(mask); colormap gray;
                
                %           pause(1);
              end
            end
            
            
            if params.doCalcBorderAngle
              %now compute the structure tensor of the local variance in phase
              minY=max(1,round(borders.borderPosY{maskfileIdx}(obj))-params.borderEvalNghSize);
              maxY=min(size(features,1),round(borders.borderPosY{maskfileIdx}(obj))+params.borderEvalNghSize);
              minX=max(1,round(borders.borderPosX{maskfileIdx}(obj))-params.borderEvalNghSize);
              maxX=min(size(features,2),round(borders.borderPosX{maskfileIdx}(obj))+params.borderEvalNghSize);
              extractids={minY:maxY,minX:maxX,':'};
              locAct=features(extractids{:});
              if ~params.borderAngleUseActivityWeighting
                locAct = ones(size(locAct));
              end
              sumLocAct = sum(locAct,3);
              if any(sumLocAct(:)==0)
                locAct = bsxfun(@plus,locAct,(sumLocAct==0)*1e-10); %add small numbers at points where all activity==0
                sumLocAct = sum(locAct,3);
              end
              sumLocPhase = sum(locAct.*exp(1i*phase(extractids{:})),3);
              
              locKernel = [0 1 0; 1 1 1; 0 1 0];
              sumLocPhase = conv2(sumLocPhase,locKernel,'same');
              sumLocAct = conv2(sumLocAct,locKernel,'same');
              locCircVar = 1 - abs(sumLocPhase./sumLocAct);
              
              diffKernel=[-0.5 0 0.5];
              Ix = conv2(locCircVar,diffKernel,'same');
              Iy = conv2(locCircVar,diffKernel','same');
              
              Ixx = Ix.^2;
              Ixy = Ix.*Iy;
              Iyy = Iy.^2;
              
              sigma=params.borderAngleSigma;
              smoothkernel=SegmentationEval.gaussian2d(size(Ix,1),size(Ix,2),sigma);
              
              meanIxx=sum(smoothkernel(:).*Ixx(:));
              meanIxy=sum(smoothkernel(:).*Ixy(:));
              meanIyy=sum(smoothkernel(:).*Iyy(:));
              
              [V,D] = eig([meanIxx meanIxy; meanIxy meanIyy]);
              if D(1,1)<D(2,2)
                structAngle=atan2(V(2,1),V(1,1));
              else
                structAngle=atan2(V(2,2),V(1,2));
              end
              structAngle=mod(structAngle,pi);
              labeledAngle=mod(borders.borderAngle{maskfileIdx}(obj),pi);
              angleError=mod(structAngle-labeledAngle,pi);
              if angleError>pi/2
                angleError=angleError-pi;
              end
              angleError=abs(angleError);
              outputBorders{simIdx,maskfileIdx}{obj}.labeledAngle = labeledAngle;
              outputBorders{simIdx,maskfileIdx}{obj}.structAngle = structAngle;
              outputBorders{simIdx,maskfileIdx}{obj}.angleError = angleError;
              if maskfileIdx==simIdx
                tmpAngleErrorMatching = tmpAngleErrorMatching+angleError;
                tmpAngleErrorMatchingCounter = tmpAngleErrorMatchingCounter+1;
              else
                tmpAngleErrorNonMatching = tmpAngleErrorNonMatching+angleError;
                tmpAngleErrorNonMatchingCounter = tmpAngleErrorNonMatchingCounter+1;
              end
              
              if params.verboseLevel>2
                disp(['structAngle:' num2str(180*structAngle/pi)]);
                disp(['labeledAngle:' num2str(180*labeledAngle/pi)]);
                disp(['angleError:' num2str(180*angleError/pi)]);
                
                figure(4)
                imagesc(locCircVar); colormap gray;
                pause(1);
              end
              
            end
          end
          
        end
        
      end
      
      
      %% Now compare matching and nonmatching:
      emptyMaskIds = cellfun(@(x) isempty(x), segmentSize);
      disp(['sum(emptyMaskIds(:)) = ' num2str(sum(emptyMaskIds(:)))]);
      
      if sum(any(~cellfun(@isempty,output)))
        [allMatching, allNonmatching] = SegmentationEval.parseResults(output);
        fnames=fieldnames(allMatching);
        for fname=1:length(fnames)
          results.matching.(fnames{fname}) = cell2mat({allMatching.(fnames{fname})});
          results.nonmatching.(fnames{fname}) = cell2mat({allNonmatching.(fnames{fname})});
          
          results.meanMatching.(fnames{fname}) = mean(results.matching.(fnames{fname}));
          results.meanNonmatching.(fnames{fname}) = mean(results.nonmatching.(fnames{fname})(:));
          
          results.meanMatchingSEM.(fnames{fname}) = std(results.matching.(fnames{fname}))/sqrt(length(results.matching.(fnames{fname})));
          results.meanNonmatchingSESegmentationEvalM.(fnames{fname}) = std(results.nonmatching.(fnames{fname})(:))/sqrt(length(results.nonmatching.(fnames{fname})(:)));
        end
      end
      
      if sum(any(~cellfun(@isempty,outputBorders)))
        [allMatchingBorders, allNonmatchingBorders] = SegmentationEval.parseResults(outputBorders);
        fnames=fieldnames(allMatchingBorders);
        for fname=1:length(fnames)
          results.matching.(fnames{fname}) = cell2mat({allMatchingBorders.(fnames{fname})});
          results.nonmatching.(fnames{fname}) = cell2mat({allNonmatchingBorders.(fnames{fname})});
          
          results.meanMatching.(fnames{fname}) = mean(results.matching.(fnames{fname}));
          results.meanNonmatching.(fnames{fname}) = mean(results.nonmatching.(fnames{fname})(:));
          
          results.meanMatchingSEM.(fnames{fname}) = std(results.matching.(fnames{fname}))/sqrt(length(results.matching.(fnames{fname})));
          results.meanNonmatchingSEM.(fnames{fname}) = std(results.nonmatching.(fnames{fname})(:))/sqrt(length(results.nonmatching.(fnames{fname})(:)));
        end
      end
      
      %% save results:
      savepath = fullfile(this.workpath,params.outSegEvalFolder);
      if this.numJobs > 1
        savepath = fullfile(savepath,num2str(this.currJobid));
      end
      mkdir(savepath);
      save(fullfile(savepath,'segEval.mat'),'results','segmentSize')
      
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function [ output ] = calcSyncValues( this, activity, phase, segMask, nghMask )
      %CALCSYNCHRONY Calculate synchrony within population A or between
      %population A and B
      %   Detailed explanation goes here
      
      numPairs = this.params.SegmentationEval.numPairs;
      
      segMask=repmat(segMask,[1 1 size(phase,3)]);
      nghMask=repmat(nghMask,[1 1 size(phase,3)]);
      
      segMask(activity==0) = false;
      nghMask(activity==0) = false;
      
      fSeg=activity(segMask);
      fNgh=activity(nghMask);
      pSeg=phase(segMask);
      pNgh=phase(nghMask);
      vSeg=exp(1i*pSeg);
      vNgh=exp(1i*pNgh);
      
      if ~isempty(numPairs) && numPairs>0
        
        %% compute pairwiseSync of Segment to Neighborhood:
        idSeg = randi(numel(fSeg),[numPairs 1]);
        idNgh = randi(numel(fNgh),[numPairs 1]);
        output.meanPairwiseSyncSegToNgh = sum(fSeg(idSeg).*fNgh(idNgh).*abs( vSeg(idSeg)+vNgh(idNgh) )) / ( 2*sum(fSeg(idSeg).*fNgh(idNgh)) );
        output.meanPairwiseLinearSyncSegToNgh = sum(fSeg(idSeg).*fNgh(idNgh).*(pi-abs(mod(pSeg(idSeg)-pNgh(idNgh),2*pi)-pi)) ) / ( sum(fSeg(idSeg).*fNgh(idNgh)) );
        
        %% compute pairwiseSync within Segment:
        idSeg1 = randi(numel(fSeg),[numPairs 1]);
        idSeg2 = randi(numel(fSeg),[numPairs 1]);
        output.meanPairwiseSyncSeg = sum(fSeg(idSeg1).*fSeg(idSeg2).*abs( vSeg(idSeg1)+vSeg(idSeg2) )) / ( 2*sum(fSeg(idSeg1).*fSeg(idSeg2)) );
        output.meanPairwiseLinearSyncSeg = sum(fSeg(idSeg1).*fSeg(idSeg2).*(pi-abs(mod(pSeg(idSeg1)-pSeg(idSeg2),2*pi)-pi)) ) / ( sum(fSeg(idSeg1).*fSeg(idSeg2)) );
        
        %% compute pairwiseSync within (Segment+Neighborhood):
        segOrNgh = logical(randi(2,[2*numPairs 1])-1);
        idsSeg = randi(numel(fSeg),[sum( segOrNgh) 1]);
        idsNgh = randi(numel(fNgh),[sum(~segOrNgh) 1]);
        v = zeros([numPairs 2]);
        v( segOrNgh) = vSeg(idsSeg);
        v(~segOrNgh) = vNgh(idsNgh);
        p = zeros([numPairs 2]);
        p( segOrNgh) = pSeg(idsSeg);
        p(~segOrNgh) = pNgh(idsNgh);
        f = zeros([numPairs 2]);
        f( segOrNgh) = fSeg(idsSeg);
        f(~segOrNgh) = fNgh(idsNgh);
        weight = f(:,1).*f(:,2);
        totalWeight = sum(weight);
        output.meanPairwiseSyncSegUnionNgh = sum(weight.*abs( v(:,1)+v(:,2) )) / ( 2*totalWeight );
        
        %% Nonlinear Phase distance measure:
        ids = ( abs( v(:,1)+v(:,2) )/2 > output.meanPairwiseSyncSegToNgh );
        output.probPairwiseSegToNghAsync = sum(weight(ids)) / totalWeight;
        
        ids = ( abs( v(:,1)+v(:,2) )/2 < output.meanPairwiseSyncSeg );
        output.probPairwiseSegSync = sum(weight(ids)) / totalWeight;
        
        %% Linear Phase distance measure:
        ids = (  pi-abs(mod(p(:,1)-p(:,2),2*pi)-pi) < output.meanPairwiseLinearSyncSegToNgh );
        output.probPairwiseSegToNghAsyncLinear = sum(weight(ids)) / totalWeight;
        
        ids = ( pi-abs(mod(p(:,1)-p(:,2),2*pi)-pi) > output.meanPairwiseLinearSyncSeg );
        output.probPairwiseSegSyncLinear = sum(weight(ids)) / totalWeight;
      end
      
      %%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%% Sync using Subpopulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      numSubpop = this.params.SegmentationEval.numSubpop;
      numSamplesInSubpop = this.params.SegmentationEval.numSamplesInSubpop;
      if ~isempty(numSubpop) && ~isempty(numSamplesInSubpop) && numSubpop>0 && numSamplesInSubpop>0
        
        if numSamplesInSubpop>min(numel(fSeg),numel(fNgh))
          numSamplesInSubpop = min(numel(fSeg),numel(fNgh));
        end
        
        %% Calc Sync within Segment using Subpopulations:
        syncSegSubpop = zeros(1,numSubpop);
        for k=1:numSubpop
          %     idSeg = randi(numel(fSeg),[numSamplesInSubpop 1]);
          idSeg = randsubset(numel(fSeg),numSamplesInSubpop);
          syncSegSubpop(k) = abs(sum( fSeg(idSeg).*vSeg(idSeg) )) / sum( fSeg(idSeg) ) ;
        end
        output.meanSyncSeg = mean(syncSegSubpop);
        
        %% Calc Sync within (Segment+Neighborhood) using Subpopulations:
        syncSegUnionNgh = zeros(1,numSubpop);
        for k=1:numSubpop
          segOrNgh = logical(randi(2,[numSamplesInSubpop 1])-1);
          %     idsSeg = randi(numel(fSeg),[sum( segOrNgh) 1]);
          %     idsNgh = randi(numel(fNgh),[sum(~segOrNgh) 1]);
          idsSeg = randsubset(numel(fSeg),sum( segOrNgh));
          idsNgh = randsubset(numel(fNgh),sum(~segOrNgh));
          v = zeros([numSamplesInSubpop 1]);
          v( segOrNgh) = vSeg(idsSeg);
          v(~segOrNgh) = vNgh(idsNgh);
          f = zeros([numSamplesInSubpop 1]);
          f( segOrNgh) = fSeg(idsSeg);
          f(~segOrNgh) = fNgh(idsNgh);
          syncSegUnionNgh(k) = abs(sum( f .* v )) / sum( f ) ;
        end
        output.meanSyncSegUnionNgh = mean(syncSegUnionNgh);
        
        ids = ( syncSegUnionNgh < output.meanSyncSeg );
        output.probSegSync = sum(ids) / length(syncSegUnionNgh);
        
      end
      
      output.ph = abs(sum( fSeg.*vSeg) ) / sum(fSeg);
      
    end
%     
%     function [ output ] = calcProbSegSync( this, activity, phase, segMask, nghMask )
%       %CALCSYNCHRONY Calculate synchrony within population A or between
%       %population A and B
%       %   Detailed explanation goes here
%       
%       numSubpop = this.params.SegmentationEval.numSubpop;
%       
%       segMask=repmat(segMask,[1 1 size(phase,3)]);
%       nghMask=repmat(nghMask,[1 1 size(phase,3)]);
%       
%       fSeg=activity(segMask);
%       fNgh=activity(nghMask);
%       pSeg=phase(segMask);
%       pNgh=phase(nghMask);
%       vSeg=exp(1i*pSeg);
%       vNgh=exp(1i*pNgh);
%       
%       %%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %%%%%%%%%%%%%%%% Sync using Subpopulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       
%       %% Calc Sync within Segment:
%       output.ph = abs(sum( fSeg.*vSeg) ) / sum(fSeg);
%       
%       %% Calc Sync within (Segment+Neighborhood) using Subpopulations:
%       if numSubpop>0
%         numSamplesInSubpop = numel(fSeg); %% this is a bug?!
%         syncSegUnionNgh = zeros(1,numSubpop);
%         for k=1:numSubpop
%           segOrNgh = logical(randi(2,[numSamplesInSubpop 1])-1);
%           idsSeg = randsubset(numel(fSeg),sum( segOrNgh));
%           idsNgh = randsubset(numel(fNgh),sum(~segOrNgh));
%           v = zeros([numSamplesInSubpop 1]);
%           v( segOrNgh) = vSeg(idsSeg);
%           v(~segOrNgh) = vNgh(idsNgh);
%           f = zeros([numSamplesInSubpop 1]);
%           f( segOrNgh) = fSeg(idsSeg);
%           f(~segOrNgh) = fNgh(idsNgh);
%           syncSegUnionNgh(k) = abs(sum( f .* v )) / sum( f ) ;
%         end
%         output.meanSyncSegUnionNgh = mean(syncSegUnionNgh);
%         
%         ids = ( syncSegUnionNgh < output.ph );
%         output.probSegSync = sum(ids) / length(syncSegUnionNgh);
%       end
%       
%     end

    function [data] = getPlotdata(segEvalObj, segSize)
      params = segEvalObj.params.SegmentationEval;
      combParallel = segEvalObj.params.Gridjob.combParallel;
      
      if nargin<2 || isempty(segSize)
        segSize = [];
      end
      if length(segSize)<3
        noSegsizeAxis = true;
      else
        noSegsizeAxis = false;
      end
      
      %% load results:
      savepath = fullfile(segEvalObj.workpath,params.outSegEvalFolder);
      if segEvalObj.numJobs > 1
        results = load(fullfile(savepath,num2str(1),'segEval.mat'),'results','segmentSize');
      else
        results = load(fullfile(savepath,'segEval.mat'),'results','segmentSize');
      end
      results = results.results;
      
      if isfield(results.matching,'meanSyncSeg')
        results.matching.segIndexNew = results.matching.meanSyncSeg-results.matching.meanSyncSegUnionNgh;
        results.nonmatching.segIndexNew = results.nonmatching.meanSyncSeg-results.nonmatching.meanSyncSegUnionNgh;
        results.meanMatching.segIndexNew = results.meanMatching.meanSyncSeg-results.meanMatching.meanSyncSegUnionNgh;
        results.meanNonmatching.segIndexNew = results.meanNonmatching.meanSyncSeg-results.meanNonmatching.meanSyncSegUnionNgh;
      end
      
      fnames=fieldnames(results.meanMatching);
      
      if segEvalObj.numJobs==1
        axisLabels{1} = 'varparams';
        axisTicks = cell(segEvalObj.numJobs,1);
        axisDims = 1;
      else
        numVarParams = length(segEvalObj.variableParams);
        if numVarParams>1 && combParallel
          axisLabels{1} = 'varparams';
          axisTicks{1} = cell(segEvalObj.numJobs,1);
          axisDims = segEvalObj.numJobs;
          for k=1:segEvalObj.numJobs
            for j=1:numVarParams
              axisTicks{1}{k} = [axisTicks{1}{k} segEvalObj.variableParams{j}{2} '=' num2str(segEvalObj.paramComb(j,k)) ' '];
            end
            axisTicks{1}{k} = axisTicks{1}{k}(1:end-1);
          end
        else
          axisLabels = cell(numVarParams,1);
          axisTicks = cell(numVarParams,1);
          axisDims = zeros(numVarParams,1);
          for j=1:numVarParams
            axisLabels{j} = segEvalObj.variableParams{j}{2};
            axisTicks{j} = segEvalObj.params.(segEvalObj.variableParams{j}{1}).(segEvalObj.variableParams{j}{2});
            axisDims(j) = length(axisTicks{j});
          end
        end
      end
      
      if ~noSegsizeAxis
        axisLabels{end+1} = 'Segment Size [ #pixels ]';
        axisTicks{end+1} = num2cell(segSize(1:end-1)+diff(segSize)/2);
        axisDims(end+1) = length(segSize)-1;
      end
      
      if length(axisDims)==1
        axisDims(end+1)=1;
      end
      
      for fname=1:length(fnames)
        meanMatching = zeros(segEvalObj.numJobs,length(segSize)-1);
        stdMatching = zeros(segEvalObj.numJobs,length(segSize)-1);
        meanNonmatching = zeros(segEvalObj.numJobs,length(segSize)-1);
        stdNonmatching = zeros(segEvalObj.numJobs,length(segSize)-1);
        meanMatchingSEM = zeros(segEvalObj.numJobs,length(segSize)-1);
        meanNonmatchingSEM = zeros(segEvalObj.numJobs,length(segSize)-1);
        pairedDiff = zeros(segEvalObj.numJobs,length(segSize)-1);
        pairedDiffCiLow = zeros(segEvalObj.numJobs,length(segSize)-1);
        pairedDiffCiUp = zeros(segEvalObj.numJobs,length(segSize)-1);
        numSegs = zeros(segEvalObj.numJobs,length(segSize)-1);
        
        breakfor = false;
        for paramid=1:segEvalObj.numJobs
          if segEvalObj.numJobs > 1
            results = load(fullfile(savepath,num2str(paramid),'segEval.mat'),'results','segmentSize');
          else
            results = load(fullfile(savepath,'segEval.mat'),'results','segmentSize');
          end
          segmentSize = cell2mat(results.segmentSize);
          results = results.results;
          
          if isfield(results.matching,'meanSyncSeg')
            results.matching.segIndexNew = results.matching.meanSyncSeg-results.matching.meanSyncSegUnionNgh;
            results.nonmatching.segIndexNew = results.nonmatching.meanSyncSeg-results.nonmatching.meanSyncSegUnionNgh;
            results.meanMatching.segIndexNew = results.meanMatching.meanSyncSeg-results.meanMatching.meanSyncSegUnionNgh;
            results.meanNonmatching.segIndexNew = results.meanNonmatching.meanSyncSeg-results.meanNonmatching.meanSyncSegUnionNgh;
          end
          
          if ~isempty(segSize)
            [numSegsTmp,bin] = histc(segmentSize,segSize);
            numSegs(paramid,:)=numSegsTmp(1:end-1);
            bin(bin==length(segSize))=0;
          else
            bin=ones(size(results.matching.(fnames{fname})));
          end
          
          if length(bin)~=length(results.matching.(fnames{fname}))
            bin=ones(size(results.matching.(fnames{fname})));
%             breakfor = true;
%             break;
          end
          
          for j=max(min(bin),1):max(bin)
            segIds = bin==j;
            matchingThisSize=results.matching.(fnames{fname})(segIds);
            tmp = results.nonmatching.(fnames{fname});
            tmp = reshape(tmp,[length(results.matching.(fnames{fname})) numel(tmp)/length(results.matching.(fnames{fname}))])';
            nonmatchingThisSize=tmp(:,segIds);
            
            %not paired:
            meanMatching(paramid,j) = mean(matchingThisSize);
            stdMatching(paramid,j) = std(matchingThisSize);
            meanNonmatching(paramid,j) = mean(nonmatchingThisSize(:));
            stdNonmatching(paramid,j) = std(nonmatchingThisSize(:));
            meanMatchingSEM(paramid,j) = std(matchingThisSize)/sqrt(length(matchingThisSize));
            meanNonmatchingSEM(paramid,j) = std(nonmatchingThisSize(:))/sqrt(length(nonmatchingThisSize(:)));
            
            %paired difference:
            matchingThisSize = repmat(matchingThisSize,[size(nonmatchingThisSize,1) 1]);
            pairedDiff(paramid,j) = mean(matchingThisSize(:)-nonmatchingThisSize(:));
            [~,~,ci]=ttest(matchingThisSize(:),nonmatchingThisSize(:));
            pairedDiffCiLow(paramid,j)=ci(1);
            pairedDiffCiUp(paramid,j)=ci(2);
          end
        end
        
%         if breakfor
%           continue;
%         end
        
%         if ~combParallel
%           if noSegsizeAxis
%             %check how many vars per param:
%             newSize = ones(1,numVarParams);
%             for j=1:numVarParams
%               newSize(j) = length(segEvalObj.params.(segEvalObj.variableParams{j}{1}).(segEvalObj.variableParams{j}{2}));
%             end
% 
%             pairedDiff = permute(reshape(pairedDiff,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
%             pairedDiffCiLow = permute(reshape(pairedDiffCiLow,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
%             pairedDiffCiUp = permute(reshape(pairedDiffCiUp,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
%             meanMatching = permute(reshape(meanMatching,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
%             meanNonmatching = permute(reshape(meanNonmatching,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
%             meanMatchingSEM = permute(reshape(meanMatchingSEM,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
%             meanNonmatchingSEM = permute(reshape(meanNonmatchingSEM,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
%             
%             pairedDiff = reshape(pairedDiff,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
%             pairedDiffCiLow = reshape(pairedDiffCiLow,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
%             pairedDiffCiUp = reshape(pairedDiffCiUp,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
%             meanMatching = reshape(meanMatching,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
%             meanNonmatching = reshape(meanNonmatching,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
%             meanMatchingSEM = reshape(meanMatchingSEM,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
%             meanNonmatchingSEM = reshape(meanNonmatchingSEM,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
%           end
%         end
        
        data.values.(fnames{fname}).pairedDiff=reshape(pairedDiff,axisDims);
        data.values.(fnames{fname}).pairedDiffCiLow=reshape(pairedDiffCiLow,axisDims);
        data.values.(fnames{fname}).pairedDiffCiUp=reshape(pairedDiffCiUp,axisDims);
        data.values.(fnames{fname}).meanMatching=reshape(meanMatching,axisDims);
        data.values.(fnames{fname}).stdMatching=reshape(stdMatching,axisDims);
        data.values.(fnames{fname}).meanNonmatching=reshape(meanNonmatching,axisDims);
        data.values.(fnames{fname}).stdNonmatching=reshape(stdNonmatching,axisDims);
        data.values.(fnames{fname}).meanMatchingSEM=reshape(meanMatchingSEM,axisDims);
        data.values.(fnames{fname}).meanNonmatchingSEM=reshape(meanNonmatchingSEM,axisDims);
        if ~isempty(numSegs)
          data.values.(fnames{fname}).numSegs=reshape(numSegs,axisDims);
        end
        
      end
      
      data.axisLabels = axisLabels;
      data.axisTicks = axisTicks;
      data.axisDims = axisDims;
      
    end
    
    function plot(segEvalObj, segSize, plotPairedDiff, noVarparamAxis)
      
      params = segEvalObj.params.SegmentationEval;
      
      combParallel = segEvalObj.params.Gridjob.combParallel;
      
      if nargin<2 || isempty(segSize)
        segSize = [0 Inf];
      end
      
      if nargin<3 || isempty(plotPairedDiff)
        plotPairedDiff = true;
      end
      
      if nargin<4 || isempty(noVarparamAxis)
        noVarparamAxis = false;
      end
      
      if length(segSize)<3
        noSegsizeAxis = true;
      else
        noSegsizeAxis = false;
      end
      
      %% load results:
      savepath = fullfile(segEvalObj.workpath,params.outSegEvalFolder);
      if segEvalObj.numJobs > 1
        results = load(fullfile(savepath,num2str(1),'segEval.mat'),'results','segmentSize');
      else
        results = load(fullfile(savepath,'segEval.mat'),'results','segmentSize');
      end
      results = results.results;
      
      fnames=fieldnames(results.meanMatching);
      
      [Selection,ok] = listdlg('ListString',fnames,'PromptString','Select fields to plot');
      fnames = fnames(Selection);
      
      
      numEvals = size(segEvalObj.paramComb,2);
      numVarParams = length(segEvalObj.variableParams);
      
      if numVarParams<=1
        combParallel = true;
      end
        
      varparamIds = 1:size( segEvalObj.paramComb,2);
      
      [data] = getPlotdata(segEvalObj, segSize);
      
      for j=1:length(fnames)
        
        if noVarparamAxis
          xaxisValues = segSize(1:end-1)+diff(segSize)/2;
          xaxisString = 'Segment Size [ #pixels ]';
          for k=1:size( segEvalObj.paramComb,1)
            tmpParamComb(k,:) = strcat(segEvalObj.variableParams{k}{2},'=',cellfun(@(x) num2str(x),segEvalObj.paramComb(k,:),'UniformOutput',false),',');
          end
          tmpParamComb = cellfun(@(x) x(1:end-1),tmpParamComb(end,:),'UniformOutput',false);
          for k=1:size( segEvalObj.paramComb,2)
            legendentries{k} = strcat(tmpParamComb{:,k});
          end
          legendentries = legendentries(varparamIds);
          
          if plotPairedDiff
            figure;
            errorbar(repmat(xaxisValues,[size(pairedDiff,1) 1])',pairedDiff',pairedDiff'-pairedDiffCiLow',-pairedDiff'+pairedDiffCiUp');
            set(gca,'XScale','log')
            set(gca,'fontsize',14,'linewidth',1);
            ylabel(['Paired Diff (' fnames{fname} ')'],'fontsize',14)
            xlabel(xaxisString,'fontsize',14)
            legend(legendentries(:));
          else
            figure; clf;
            hold on;
            errorbar(repmat(xaxisValues,[size(meanMatching,1) 1])',meanMatching',meanMatchingSEM');
            errorbar(repmat(xaxisValues,[size(meanNonmatching,1) 1])',meanNonmatching',meanNonmatchingSEM','-.');
            hold off;
            set(gca,'XScale','log')
            set(gca,'fontsize',14,'linewidth',1);
            xlabel(xaxisString,'fontsize',14)
            ylabel(fnames{fname},'fontsize',14)
            legendentries=strcat(repmat(cellfun(@(x) num2str(x),legendentries','UniformOutput',false),[1 2]),repmat({' '},[size(legendentries,2) 2]),repmat({'Matching','Nonmatching'},[size(legendentries,2) 1]));
            legend(legendentries(:));
          end
          set(gca,'XScale','log')
        elseif noSegsizeAxis
          if combParallel
            if numVarParams==0
              xaxisValues = 1;
              xaxisString = 'varparamIds';
            elseif numVarParams==1
              xaxisValues = cell2mat(segEvalObj.paramComb(1,varparamIds));
              xaxisString = segEvalObj.variableParams{1}{2};
            else
              xaxisValues = paramsString;
              xaxisString = 'varparamIds';
            end
            legendentries = {''};
          else
            xaxisValues = cell2mat(segEvalObj.params.(segEvalObj.variableParams{paramXaxisId}{1}).(segEvalObj.variableParams{paramXaxisId}{2}));
            xaxisString = segEvalObj.variableParams{paramXaxisId}{2};
            yaxisParams = cellfun(@(x) x{2}, segEvalObj.variableParams([1:paramXaxisId-1 paramXaxisId+1:end]),'UniformOutput',false);
            for k=1:size( segEvalObj.paramComb,1)
              tmpParamComb(k,:) = strcat(segEvalObj.variableParams{k}{2},'=',cellfun(@(x) num2str(x),segEvalObj.paramComb(k,:),'UniformOutput',false),',');
            end
            tmpParamComb(paramXaxisId,:) = [];
            tmpParamComb = cellfun(@(x) x(1:end-1),tmpParamComb(end,:),'UniformOutput',false);
            for k=1:size( segEvalObj.paramComb,2)
              tmpParamComb2{k} = strcat(tmpParamComb{:,k});
            end
            tmpParamComb2 = permute(reshape(tmpParamComb2,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
            tmpParamComb2 = reshape(tmpParamComb2,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
            legendentries = tmpParamComb2(1,:);
          end
          
          if plotPairedDiff
            
            figure;
            errorbar(repmat(xaxisValues,[size(pairedDiff,2) 1])',pairedDiff,pairedDiff-pairedDiffCiLow,-pairedDiff+pairedDiffCiUp);
            set(gca,'fontsize',14,'linewidth',1);
            xlabel(xaxisString,'fontsize',14)
            ylabel(['Paired Diff (' fnames{fname} ')'],'fontsize',14)
            
            if ~isempty(legendentries) && ~sum(strcmp(legendentries{1},''))
              legend(legendentries(:));
            end

          else
            
            figure; clf;
            hold on;
            errorbar(repmat(xaxisValues,[size(meanMatching,2) 1])',meanMatching,meanMatchingSEM);
            errorbar(repmat(xaxisValues,[size(meanNonmatching,2) 1])',meanNonmatching,meanNonmatchingSEM,'-.');
            hold off;
            set(gca,'fontsize',14,'linewidth',1);
            xlabel(xaxisString,'fontsize',14)
            ylabel(fnames{fname},'fontsize',14)
            
            legendentries=strcat(repmat(cellfun(@(x) num2str(x),legendentries','UniformOutput',false),[1 2]),repmat({' '},[size(legendentries,2) 2]),repmat({'Matching','Nonmatching'},[size(legendentries,2) 1]));
            legend(legendentries(:));
          
          end
        else
          yaxisValues = segSize(1:end-1)+diff(segSize)/2;
          
          if combParallel
            if numVarParams==0
              xaxisValues = 1;
              xaxisString = 'varparamIds';
            elseif numVarParams==1
              xaxisValues = cell2mat(segEvalObj.paramComb(1,varparamIds));
              xaxisString = segEvalObj.variableParams{1}{2};
            else
              xaxisValues = varparamIds;
              xaxisString = 'varparamIds';
            end
            legendentries = {''};
          else
            xaxisValues = cell2mat(segEvalObj.params.(segEvalObj.variableParams{paramXaxisId}{1}).(segEvalObj.variableParams{paramXaxisId}{2}));
            xaxisString = segEvalObj.variableParams{paramXaxisId}{2};
            yaxisParams = cellfun(@(x) x{2}, segEvalObj.variableParams([1:paramXaxisId-1 paramXaxisId+1:end]),'UniformOutput',false);
            for k=1:size( segEvalObj.paramComb,1)
              tmpParamComb(k,:) = strcat(segEvalObj.variableParams{k}{2},'=',cellfun(@(x) num2str(x),segEvalObj.paramComb(k,:),'UniformOutput',false),',');
            end
            tmpParamComb(paramXaxisId,:) = [];
            tmpParamComb = cellfun(@(x) x(1:end-1),tmpParamComb(end,:),'UniformOutput',false);
            for k=1:size( segEvalObj.paramComb,2)
              tmpParamComb2{k} = strcat(tmpParamComb{:,k});
            end
            tmpParamComb2 = permute(reshape(tmpParamComb2,newSize),[paramXaxisId 1:paramXaxisId-1 paramXaxisId+1:length(newSize)]);
            tmpParamComb2 = reshape(tmpParamComb2,[size(pairedDiff,1) numel(pairedDiff)/size(pairedDiff,1)]);
            legendentries = tmpParamComb2(1,:);
          end
          
          
          if plotPairedDiff
            figure;
            surf(yaxisValues,xaxisValues,pairedDiff)
            set(gca,'XScale','log')
            set(gca,'fontsize',14,'linewidth',1);
            xlabel('Segment Size  [ #pixels ]','fontsize',14)
            ylabel(xaxisString,'fontsize',14)
            zlabel(['Paired Diff (' fnames{fname} ')'],'fontsize',14)
          else
            figure;
            surf(yaxisValues,xaxisValues,meanMatching)
            set(gca,'XScale','log')
            set(gca,'fontsize',14,'linewidth',1);
            xlabel('Segment Size  [ #pixels ]','fontsize',14)
            ylabel(xaxisString,'fontsize',14)
            zlabel(['Matching ' fnames{fname}],'fontsize',14)
            
            figure;
            surf(yaxisValues,xaxisValues,meanNonmatching)
            set(gca,'XScale','log')
            set(gca,'fontsize',14,'linewidth',1);
            xlabel('Segment Size  [ #pixels ]','fontsize',14)
            ylabel(xaxisString,'fontsize',14)
            zlabel(['Nonmatching ' fnames{fname}],'fontsize',14)
          end
          
        end
        
%         if plotNumSeg
%           figure;
%           surf(yaxisValues,varparamIds,numSegs)
%           set(gca,'XScale','log')
%           set(gca,'fontsize',14,'linewidth',1);
%           xlabel('Segment Size  [ #pixels ]','fontsize',14)
%           ylabel('varparamIds','fontsize',14)
%           zlabel('Number of Segments','fontsize',14)
%         end
        
%         if nargin>4 && ~isempty(outFilename)
%           matlab2tikz( [outFilename '_' fnames{fname} '.tikz'] , 'width', '\figurewidth', 'relativePngPath', '', 'interpretTickLabelsAsTex', true );
%         end
      end
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this,segSizeBinEdges)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
%       if nargin>1
%         this.params.SegmentationEval.segSizeBinEdges = segSizeBinEdges;%exp(5:10);
%       end
%       
%       dataPerSegSize = getPlotdata(this, this.params.SegmentationEval.segSizeBinEdges);
%       data = getPlotdata(this);
%       
%       paramComb = this.paramComb;
%       variableParams = this.variableParams;
%       
%       savepath = fullfile(this.resultpath,this.params.Gridjob.jobname);
%       mkdir(savepath);
%       save(fullfile(savepath,'segEval.mat'),'dataPerSegSize','data','paramComb','variableParams');
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
  methods (Static)
    function [allMatching, allNonmatching] = parseResults(resultOutput)
      %% calc resultMatching and meanResultNonmatching:
      objCounter=1;
      for maskFileIdx=1:size(resultOutput,1)
        for objIdx=1:length(resultOutput{maskFileIdx,maskFileIdx})
          allMatching(objCounter,1) = resultOutput{maskFileIdx,maskFileIdx}{objIdx};
          nonmatchCounter=1;
          for simIdx=1:size(resultOutput,1)
            if simIdx~=maskFileIdx && ~isempty(resultOutput{simIdx,maskFileIdx}{objIdx})
              allNonmatching(objCounter,nonmatchCounter) = resultOutput{simIdx,maskFileIdx}{objIdx};
              nonmatchCounter=nonmatchCounter+1;
            end
          end
          objCounter=objCounter+1;
        end
      end
    end
    
    function f=gaussian2d(Nx,Ny,sigma)
      % N is grid size, sigma speaks for itself
      Nx=Nx-1;
      Ny=Ny-1;
      [x y]=meshgrid(-Nx/2:Nx/2, -Ny/2:Ny/2);
      f=exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
      f=f'./sum(f(:));
    end
  end
  
end

