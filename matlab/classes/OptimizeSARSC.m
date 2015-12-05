classdef OptimizeSARSC < Gridjob
  %OptimizeSARSC Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = OptimizeSARSC(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
            
      this.params.OptimizeSARSC.constrainSym = false; % if yes SC is just symmetric, if 'rowNorm' then SC will first be row normalized and then transpose is added
      this.params.OptimizeSARSC.useFrob = false;
      this.params.OptimizeSARSC.useRowRenorm = 'no';  % or 'yes' or 'allowScaling'  
      this.params.OptimizeSARSC.constrainPos = false;
      this.params.OptimizeSARSC.k = 0.65;
      this.params.OptimizeSARSC.lambda = 0;
      this.params.OptimizeSARSC.gamma = 1;
      this.params.OptimizeSARSC.negPenalty = 0;
      this.params.OptimizeSARSC.maxSteps = 65536;
      this.params.OptimizeSARSC.savefolder = 'results';
      this.params.OptimizeSARSC.unittest = false;
      this.params.OptimizeSARSC.useRMSprop = true;
      this.params.OptimizeSARSC.rmsPropDecay = 0.9;
      this.params.OptimizeSARSC.initLearnRate = 1e-4;
      this.params.OptimizeSARSC.saveAtIters = [0 2.^(0:16)];
      this.params.OptimizeSARSC.breakAtCost = 0.999;
      this.params.OptimizeSARSC.noiseAmount = 0.5;
      this.params.OptimizeSARSC.use_fmincon = false;
      this.params.OptimizeSARSC.loglevel = 1;
      
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
      
      p = this.params.OptimizeSARSC;
      
      paths = dataPaths();
      
      outstr = '';
      for k=1:length(this.variableParams)
        value = this.paramComb{k,this.currJobid};
        if ~ischar(value)
          value = num2str(value);
        end
        outstr = [outstr this.variableParams{k}{2} '_' value '_'];
      end
      
      disp(outstr)
      
      load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));
      empSC = avg_ci;
      empSC(isnan(empSC)) = 0;
      empSC = empSC + empSC';
      
      if p.constrainSym
        empSC = normGraph(empSC, avg_roi_size, 'ROIprd', false, 0);
        if p.constrainSym=='rowNorm'
          empSC = bsxfun(@rdivide,empSC,sum(empSC,2));
          empSC = empSC + empSC';
        end
        empSC = 66 * empSC / sum(empSC(:));
      else
        empSC = normGraph(empSC, avg_roi_size, 'ROIprd', true, 0);
      end
      
      dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150125_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30_entmirrored.mat']);
      tmp = mean(dataEegTmp.coh_all([1:13 15 17:20],:,5:6,:,:,7),3);
      tmp = mean(tmp,1);
      empFC = squeeze(tmp);
      clear dataEegTmp;
      
      triuIds = find(triu(ones(66,66),1));
      offDiagIds = find(~eye(66,66));
      
      %% outer for loop for different noise distributions
      noise = rand(size(empSC));
      noise(logical(eye(size(empSC)))) = 0;
      if p.constrainSym
        noise(triuIds) = 0; %#ok<FNDSB>
        noise = noise + noise';
        noise = 66 * noise / sum(noise(:));
      else
        noise = bsxfun(@rdivide,noise,sum(noise,2));
      end
      
      interpCurrent = p.noiseAmount;
      initSC = (1-interpCurrent)*empSC + interpCurrent*noise;
      [finalLearnSC, dev, logLearnSC, cc] = optimizeSC(initSC, empSC, empFC, p);
      finalLearnSCoffDiag = finalLearnSC(offDiagIds);
      
      tmp = find(dev==0, 1, 'first');
      if isempty(tmp)
        numIters = length(dev);
      else
        numIters = tmp-1;
      end
      
      finalCost = dev(numIters);
      
      %% evaluate results:
      mkdir(fullfile(this.workpath, p.savefolder));
      
      save(fullfile(this.workpath, p.savefolder, ['result_job' num2str(this.currJobid) '_' outstr  '.mat']));


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
    
    %% gather all the results
    function [results] = collectResults(this)
      
      p = this.params.OptimizeSARSC;
      loadpath = fullfile(this.workpath, p.savefolder);
      
      if ~exist(loadpath,'dir')
          % check if we find the files in a relative path:
          loadpath = uigetdir(pwd);
      end
      
      if exist(fullfile(loadpath,'results.mat'),'file')
        
        results = load(fullfile(loadpath,'results.mat'));
        
      else
        
        results.learnedSC = cell(1,this.numJobs);
        results.itersSaved = cell(1,this.numJobs);
        results.dev = cell(1,this.numJobs);
        results.cc = cell(1,this.numJobs);
        results.numSaved = zeros(1,this.numJobs);

        for i=1:this.numJobs
          fname=dir(fullfile(loadpath, ['result_job' num2str(i) '_*.mat']));
          tmp=load(fullfile(loadpath, fname.name));

          %% readout SC:
          learnedSC = cell2mat(permute(tmp.logLearnSC,[2 3 1]));
          numSaved = size(learnedSC,3);
          itersSaved = tmp.p.saveAtIters(1:numSaved);

          %% append last iter if it was not saved normally:
          lastIter = find(tmp.dev==0,1,'first');
          if ~isempty(lastIter)
            learnedSC(:,:,end+1) = tmp.finalLearnSC;
            numSaved = numSaved + 1;
            itersSaved(end+1) = lastIter;
          end

          %% save deviance at the corresponding positions:
          % add 1 because we gathered the dev within the next run.
          dev = tmp.dev(min(tmp.p.saveAtIters+1,length(tmp.dev)));
          cc = tmp.cc(min(tmp.p.saveAtIters+1,length(tmp.dev)));

          results.learnedSC{i} = learnedSC;
          results.numSaved(i) = numSaved;
          results.itersSaved{i} = itersSaved;
          results.dev{i} = dev;
          results.cc{i} = cc;

        end
        
        save(fullfile(loadpath,'results.mat'), '-struct', 'results');
        
      end
      
    end
    
    %% plotting results
    function plot(this,dimsSelect,distMetId,doCreateVideo,doPlotSCs,appendAnim)
      
      if nargin<2
        dimsSelect = [];
      end
      if nargin<3
        distMetId = [];
      end
      if nargin<4
        doCreateVideo = [];
      end
      if nargin<5
        doPlotSCs = true;
      end
      if nargin<6
        appendAnim = false;
      end
      
      results = this.collectResults();
      
      %% select dimensions:
      for k=1:length(this.variableParams)
          varParams{k} = this.variableParams{k}{2};
      end
      
      if ~isempty(dimsSelect)
        dim1 = dimsSelect(1);
        dim2 = dimsSelect(2);
      else
        [dim1,v] = listdlg('PromptString','Select first dimension:','SelectionMode','single','ListString',varParams);

        dimsRem = setdiff(1:length(varParams),dim1);
        varParamsRem = varParams(dimsRem);
        [dimSel,v] = listdlg('PromptString','Select second dimension:','SelectionMode','single','ListString',varParamsRem);
        dim2 = dimsRem(dimSel);
      end
      
      disp(['selected dimensions: ' num2str(dim1) ' , ' num2str(dim2)])
      
      dim1values = this.params.(this.variableParams{dim1}{1}).(this.variableParams{dim1}{2});
      dim2values = this.params.(this.variableParams{dim2}{1}).(this.variableParams{dim2}{2});
      
      dim1valuesPerJob = cell2mat(this.paramComb(dim1,:));
      dim2valuesPerJob = cell2mat(this.paramComb(dim2,:));
      
      %% 
      if doPlotSCs
        dim1and2Values = [dim1values{1} dim2values{1}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_init_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )

        dim1and2Values = [dim1values{1} dim2values{end}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_init_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )

        dim1and2Values = [dim1values{end} dim2values{1}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_init_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )

        dim1and2Values = [dim1values{end} dim2values{end}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_init_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )

        %
        dim1and2Values = [dim1values{1} dim2values{1}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_finish_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )

        dim1and2Values = [dim1values{1} dim2values{end}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_finish_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )

        dim1and2Values = [dim1values{end} dim2values{1}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_finish_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )

        dim1and2Values = [dim1values{end} dim2values{end}];
        jobId = find(dim1valuesPerJob==dim1and2Values(1) & dim2valuesPerJob==dim1and2Values(2));
        fname = ['SC_finish_' varParams{dim1} '=' num2str(dim1and2Values(1)) ' and ' varParams{dim2} '=' num2str(dim1and2Values(2))];
        plotConnMat( results.learnedSC{jobId}(:,:,1), pwd, fname, [], true, true, true, true )
      end
      
      %% convert to diag form
      offDiagIds = find(~eye(66,66));
      vecSC = cellfun(@(x) reshape(x,[66*66 size(x,3)]), results.learnedSC, 'UniformOutput', false);
      vecSC = cellfun(@(x) x(offDiagIds,:), vecSC, 'UniformOutput', false);
      vecSC_final = cell2mat( cellfun(@(x) x(:,end), vecSC, 'UniformOutput', false) );
      vecSC_init = cell2mat( cellfun(@(x) x(:,1), vecSC, 'UniformOutput', false) );
      
      %% calc similarity matrix:
      distMet = {'correlation', 'euclidean', 'cosine', 'cityblock'};
      if isempty(distMetId)
        [distMetId,v] = listdlg('PromptString','Select distance metric:','SelectionMode','multiple','ListString',distMet);
      end
      
      %%
      for d=1:length(distMetId)
        
        currentDistMet = distMet{distMetId(d)};
        
        %% calc distances between all sims and saved iterations
        allVecSCs = cell2mat(vecSC);
        allIters = cell2mat(cellfun(@(x) 1:size(x,2), vecSC, 'UniformOutput', false));
        allDev = cell2mat(results.dev')';
        allCorr = cellfun(@(x) x.corr, cat(1,results.cc{:}))';
        allJobIds = cell2mat(cellfun(@(x) x(2)*ones(1,x(1)), num2cell([results.numSaved; 1:length(results.numSaved)],1), 'UniformOutput', false));
        distances = pdist(allVecSCs',currentDistMet);


        %%
        if exist(['mds3d_' currentDistMet '.mat'],'file')
          disp('reload mds3d...')
          Y3 = load(['mds3d_' currentDistMet '.mat']);
          Y3 = Y3.Y3;
        else
          disp('calc mds3d...')
          Y3 = mdscale(distances,3,'Criterion','sstress');
          save(['mds3d_' currentDistMet '.mat'],'Y3');
        end

        %% plot 3d mds with color as loss:
        figh = figure(3);
        clf;
        this.plot_mds(Y3, allJobIds, allCorr, allIters, dim1valuesPerJob, dim2valuesPerJob, dim1values,dim2values, Inf)

        %%
        if isempty(doCreateVideo)
          doCreateVideo = questdlg('Create Video?','Video', 'Yes','No','No');
        end

        switch doCreateVideo
          case {'No',false}
            
          case {'Yes',true}
            if ispc
              daObj=VideoWriter(['mds3d_' currentDistMet '.mp4'],'MPEG-4');
            else
              daObj=VideoWriter(['mds3d_' currentDistMet '.avi']);
            end
            daObj.FrameRate=12;
            open(daObj);

            az_init = 45;
            el_init = 30;
            az = az_init;
            el = el_init;
            view(az,el);
            axis vis3d
            for t=1:180
              az=az+2;
              view(az,el);
              pause(0.05);
              drawnow;
              writeVideo(daObj,getframe(figh));
            end
            
            if appendAnim
              xlims = get(gca,'xlim');
              ylims = get(gca,'ylim');
              zlims = get(gca,'zlim');
              clims = get(gca,'clim');
              for t=1:180
                maxSaveId = 2+ceil(t/8);
                clf;
                this.plot_mds(Y3, allJobIds, allCorr, allIters, dim1valuesPerJob, dim2valuesPerJob, dim1values, dim2values, maxSaveId)
                set(gca,'xlim',xlims);
                set(gca,'ylim',ylims);
                set(gca,'zlim',zlims);
                set(gca,'clim',clims);
                
                uicontrol(gcf,'style','text','units','normalized','pos',[0.3 0.95 0.4 0.05],'string',['iteration' num2str(this.params.OptimizeSARSC.saveAtIters(maxSaveId))],'FontSize',16)

                view(az_init,el_init);
                axis vis3d

                az=az+2;
                view(az,el);
                pause(0.05);
                drawnow;
                writeVideo(daObj,getframe(figh));
              end
            end
            

            close(daObj);

        end
      end
      
    end
    
    
    function plot_mds(this, Y3, allJobIds, allDev, allIters, dim1valuesPerJob, dim2valuesPerJob, dim1values,dim2values,maxSaveId)
      
      colormap(jet)
      hold on;
      for i1=1:length(dim1values)
        for i2=1:length(dim2values)
          
          jobId = find(dim1valuesPerJob==dim1values{i1} & dim2valuesPerJob==dim2values{i2});
          SCids = find(allJobIds==jobId);
          SCids = SCids(1:min(maxSaveId,length(SCids)));
          
          %% search color:
          my_dev = allDev(SCids);
          surface([Y3(SCids,1)';Y3(SCids,1)'],[Y3(SCids,2)';Y3(SCids,2)'],[Y3(SCids,3)';Y3(SCids,3)'],[my_dev;my_dev],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2);
          
          %plot last iters
          plot3(Y3(SCids(end),1),Y3(SCids(end),2),Y3(SCids(end),3),'ro','LineWidth',3)
        end
        
        % find all first iters of all these dim1(i1) runs:
        jobIds = find(dim1valuesPerJob==dim1values{i1});
        SCids = find( ismember(allJobIds, jobIds) & allIters==1);
        plot3(Y3(SCids,1),Y3(SCids,2),Y3(SCids,3),'linew',2,'Color','k')
        plot3(Y3(SCids,1),Y3(SCids,2),Y3(SCids,3),'ko','LineWidth',3)
        
      end
      
      plot3(Y3(1,1),Y3(1,2),Y3(1,3),'go','LineWidth',3)
      %title('3D MDS')
      colorbar
%       set(gca,'clim',[0 1])
      axis equal;
      box on
      grid on
          
      az = 45;
      el = 30;
      view(az,el);
      
      set(gca,'FontSize',16')
      set(gca,'xticklabel',[])
      set(gca,'yticklabel',[])
      set(gca,'zticklabel',[])
      
      xticks = get(gca,'xtick');
      yticks = get(gca,'ytick');
      zticks = get(gca,'ztick');
      
      diff_ticks = min(min( xticks(2)-xticks(1), yticks(2)-yticks(1)), zticks(2)-zticks(1));
      
      set(gca,'xtick', -diff_ticks*10:diff_ticks:10*diff_ticks);
      set(gca,'ytick', -diff_ticks*10:diff_ticks:10*diff_ticks);
      set(gca,'ztick', -diff_ticks*10:diff_ticks:10*diff_ticks);
    end
    
  end
  
end

