classdef ConnectomeEnvelopeReduce < Gridjob
  %ConnectomeEnvelopeReduce Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ConnectomeEnvelopeReduce(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ConnectomeEnvelopeReduce.ConnectomeSimJobName = 'ConnectomeSim';
      this.params.ConnectomeEnvelopeReduce.ConnectomeSimOut = 'ConnectomeSim';
      this.params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeJobName = 'ConnectomeEnvelope';
      this.params.ConnectomeEnvelopeReduce.ConnectomeEnvelopeOut = 'connEnv';
      this.params.ConnectomeEnvelopeReduce.onlyFCsim = false; % do not use envelopes but the directly simulated FC
      
      this.params.ConnectomeEnvelopeReduce.eegDatabase = 1; 
      % 1=icoh_all_0_20140715.mat initial lcmv
      % 2=icoh_all_mne2_20140730.mat uses mni
      % 3=icoh_all_lcmvhilbert_20140804.mat includes cohy uses lcmv and hilbert
      % 4=icoh_all_lcmvhilbertrest_20140807.mat only real RS
      
      this.params.ConnectomeEnvelopeReduce.eeg = struct();
%       this.params.ConnectomeEnvelopeReduce.eeg.subj.Ids = [1:4 7:10];
%       this.params.ConnectomeEnvelopeReduce.eeg.subj.Avg = true;
%       this.params.ConnectomeEnvelopeReduce.eeg.day.Ids = [];
%       this.params.ConnectomeEnvelopeReduce.eeg.day.Avg = true;
%       this.params.ConnectomeEnvelopeReduce.eeg.prepost.Ids = 1;
%       this.params.ConnectomeEnvelopeReduce.eeg.prepost.Avg = false;
%       this.params.ConnectomeEnvelopeReduce.eeg.task.Ids = [];
%       this.params.ConnectomeEnvelopeReduce.eeg.task.Avg = true;

      this.params.ConnectomeEnvelopeReduce.sim = struct();
%       this.params.ConnectomeEnvelopeReduce.sim.Dim1.Ids = [];
%       this.params.ConnectomeEnvelopeReduce.sim.Dim1.Avg = false;
%       this.params.ConnectomeEnvelopeReduce.sim.Dim2.Ids = [];
%       this.params.ConnectomeEnvelopeReduce.sim.Dim2.Avg = false;
      
      this.params.ConnectomeEnvelopeReduce.outDirectory = 'results';
      
      this.params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no'; % 'no' or 'match' or 'nonmatch'
      
      this.params.ConnectomeEnvelopeReduce.results = struct();
%       this.params.ConnectomeEnvelopeReduce.results.Dim1.Ids = [];
%       this.params.ConnectomeEnvelopeReduce.results.Dim1.Avg = false;

      this.params.ConnectomeEnvelopeReduce.doPlot = true; 
      this.params.ConnectomeEnvelopeReduce.deleteEnvelopeWhenDone = false;
      this.params.ConnectomeEnvelopeReduce.deleteSimWhenDone = false;
      this.params.ConnectomeEnvelopeReduce.permutePlotDims = [];
      this.params.ConnectomeEnvelopeReduce.plotPermTests = true;

      this.params.ConnectomeEnvelopeReduce.reloadConnFC = false;
      this.params.ConnectomeEnvelopeReduce.reloadCompareSimExp = false;
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    %% Start: is executed before all individual parameter jobs are started
    function startJob(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
      
      
      %%%% END EDIT HERE:                                        %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algori
      
      p = this.params.ConnectomeEnvelopeReduce;
      %% for backwards compatibility:
      if isfield(p,'eegSubjIds') && ~isfield(p,'eeg') && ~isfield(p,'sim')
        this.params.ConnectomeEnvelopeReduce.eeg.subjEeg.Ids = p.eegSubjIds;
        this.params.ConnectomeEnvelopeReduce.eeg.subjEeg.Avg = true;
        this.params.ConnectomeEnvelopeReduce.eeg.day.Ids = [];
        this.params.ConnectomeEnvelopeReduce.eeg.day.Avg = true;
        this.params.ConnectomeEnvelopeReduce.eeg.prepost.Ids = 1;
        this.params.ConnectomeEnvelopeReduce.eeg.prepost.Avg = false;
        this.params.ConnectomeEnvelopeReduce.eeg.task.Ids = [];
        this.params.ConnectomeEnvelopeReduce.eeg.task.Avg = true;
        this.params.ConnectomeEnvelopeReduce.sim = struct();
      end
      
      if isfield(this.params.ConnectomeEnvelopeReduce.eeg,'subj')
        this.params.ConnectomeEnvelopeReduce.eeg.subjEeg = this.params.ConnectomeEnvelopeReduce.eeg.subj;
        this.params.ConnectomeEnvelopeReduce.eeg = rmfield(this.params.ConnectomeEnvelopeReduce.eeg,'subj');
      end
      
      if isfield(this.params.ConnectomeEnvelopeReduce.results,'subj')
        this.params.ConnectomeEnvelopeReduce.results.subjEeg = this.params.ConnectomeEnvelopeReduce.results.subj;
        this.params.ConnectomeEnvelopeReduce.results = rmfield(this.params.ConnectomeEnvelopeReduce.results,'subj');
      end
        
      savepath = fullfile(this.workpath,p.outDirectory);
      
      %% start calculations:
      
      
      if p.reloadCompareSimExp
        disp('Start loading compareSimExp...')
        compareSimExp = load(fullfile(savepath,['compareSimExp' num2str(this.currJobid)]));
      else
        
        if p.reloadConnFC
          disp('Start loading ConnFC...')
          ConnFC = load(fullfile(savepath,'ConnFC.mat'));
        else
          disp('Start calculating ConnFC...')
          ConnFC = this.gatherResults();
        end
        disp('Start calculating compareSimExp...')
        compareSimExp = this.compareWithEEG(ConnFC);
      end
      
      if p.doPlot
        disp('Start plotting...')
        this.plotResults(compareSimExp);
      end
      
      if p.deleteEnvelopeWhenDone
        disp('Start deleting Envelope directories...')
        [stat, mess, id]=rmdir([this.workpath '/temp_' p.ConnectomeEnvelopeJobName],'s');
        [stat, mess, id]=rmdir([this.workpath '/' p.ConnectomeEnvelopeOut],'s');
      end
      
      if p.deleteSimWhenDone
        disp('Start deleting Simulation directories...')
        [stat, mess, id]=rmdir([this.workpath '/temp_' p.ConnectomeSimJobName],'s');
        [stat, mess, id]=rmdir([this.workpath '/' p.ConnectomeSimOut],'s');
      end
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function ConnFC = gatherResults(this)
      
      p = this.params.ConnectomeEnvelopeReduce;
      savepath = fullfile(this.workpath,p.outDirectory);
      mkdir(savepath);
      
      varParam = load([this.workpath '/temp_' p.ConnectomeSimJobName '/jobDesc.mat'],'variableParams','paramComb','params');
      
      for j=1:length(varParam.variableParams)
        paramName{j} = varParam.variableParams{j}{2};
        paramVar{j} = varParam.params.(varParam.variableParams{j}{1}).(varParam.variableParams{j}{2});
        dims(j) = length(paramVar{j});
      end
      
      numSims = size(varParam.paramComb,2);
      
      if p.onlyFCsim
        if exist(fullfile(this.workpath, p.ConnectomeSimOut),'dir')
          %load first dataset to get dimensions:
          stats=load([this.workpath '/' p.ConnectomeSimOut '/1FC.mat']);

          numFreq = 1;
          numRois = size(stats.FCsimNoBold,1);

          plv = zeros(numRois,numRois,numFreq,numSims);
          icoh = zeros(numRois,numRois,numFreq,numSims);
          coh = zeros(numRois,numRois,numFreq,numSims);
          cohy = zeros(numRois,numRois,numFreq,numSims);
          FCsim = zeros(numRois,numRois,numFreq,numSims);
          for k=1:numSims
            stats=load([this.workpath '/' p.ConnectomeSimOut '/' num2str(k) 'FC.mat']);

            plv(:,:,1,k) = stats.FCsimNoBold;
            icoh(:,:,1,k) = stats.FCsimNoBold;
            coh(:,:,1,k) = stats.FCsimNoBold;
            cohy(:,:,1,k) = stats.FCsimNoBold;
            FCsim(:,:,1,k) = stats.FCsimNoBold;
          end
        end
      else
        %load first dataset to get dimensions:
        stats=load([this.workpath '/' p.ConnectomeEnvelopeOut '/FCsimJob1.mat']);
        numFreq = length(stats.PLVsim);
        numRois = size(stats.PLVsim{1}.PLV{1},1);

        plv = zeros(numRois,numRois,numFreq,numSims);
        icoh = zeros(numRois,numRois,numFreq,numSims);
        coh = zeros(numRois,numRois,numFreq,numSims);
        cohy = zeros(numRois,numRois,numFreq,numSims);
        FCsim = zeros(numRois,numRois,numFreq,numSims);
        meanOrderParam = zeros(1,numSims);
        stdOrderParam = zeros(1,numSims);
        for k=1:numSims
          stats=load([this.workpath '/' p.ConnectomeEnvelopeOut '/FCsimJob' num2str(k) '.mat']);
          for freq=1:numFreq
            plv(:,:,freq,k) = stats.PLVsim{freq}.PLV{1};
            icoh(:,:,freq,k) = stats.ICOHsim{freq}.icoh{1};
            coh(:,:,freq,k) = stats.ICOHsim{freq}.coherence{1};
            cohy(:,:,freq,k) = stats.ICOHsim{freq}.coherency{1};
            FCsim(:,:,freq,k) = stats.FCsim{freq}.FC{1};
          end

          if exist(fullfile(this.workpath, p.ConnectomeSimOut),'dir')
            simStats=load([this.workpath '/' p.ConnectomeSimOut '/' num2str(k) 'KuramotoStats.mat']);
            meanOrderParam(k) = simStats.meanOrderParam;
            stdOrderParam(k) = simStats.stdOrderParam;
          end
        end
      
      end
      
      plv = reshape(plv,[numRois numRois numFreq dims]);
      icoh = reshape(icoh,[numRois numRois numFreq dims]);
      coh = reshape(coh,[numRois numRois numFreq dims]);
      cohy = reshape(cohy,[numRois numRois numFreq dims]);
      FCsim = reshape(FCsim,[numRois numRois numFreq dims]);
      
      if length(dims)==1
        dims=[dims 1];
      end
      
      if exist('meanOrderParam','var')
        meanOrderParam = reshape(meanOrderParam,dims);
        stdOrderParam = reshape(stdOrderParam,dims);
        ConnFC.meanOrderParam = meanOrderParam;
        ConnFC.stdOrderParam = stdOrderParam;
      end
      
      ConnFC.plv = plv;
      ConnFC.icoh = icoh;
      ConnFC.coh = coh;
      ConnFC.cohy = cohy;
      ConnFC.FCsim = FCsim;
      
      ConnFC.paramName = paramName;
      ConnFC.paramVar = paramVar;
      ConnFC.dims = dims;
      
      save([savepath '/ConnFC.mat'],'-struct','ConnFC')
      
    end
    
    function compareSimExp = compareWithEEG(this,sim)
      p = this.params.ConnectomeEnvelopeReduce;
      savepath = fullfile(this.workpath,p.outDirectory);
      
      %% load simulation values:
      if nargin<2
        data.sim = load([savepath '/ConnFC.mat']);
      else
        data.sim = sim;
      end
      data.sim.metrics = {'icoh','coh','plv','cohy'};
      data.sim.dimName = cat(2,{'roi','roi','freq'},data.sim.paramName);
      data.sim.dimSize = size(data.sim.plv);
      data.sim.dimLabels = [{num2cell(1:66),num2cell(1:66),{5,11,25}},data.sim.paramVar];
      
      
      paths = dataPaths( );
      
      
      %% load eeg values:
      if p.eegDatabase==4
        data.eeg = load([paths.databases '/SC_Bastian/icoh_all_lcmvhilbertrest_20140807.mat']);
        data.eeg.metrics = {'icoh','coh','plv','cohy'};
        for n=1:length(data.eeg.metrics)
          % replace NaNs of subject 6 on day 2 with day 1:
          data.eeg.([data.eeg.metrics{n} '_all'])(6,2,5:6,:,:,:) = data.eeg.([data.eeg.metrics{n} '_all'])(6,1,5:6,:,:,:);
          % replace NaNs of subject 7 on day 2 post_rs with day 1 post_rs:
          data.eeg.([data.eeg.metrics{n} '_all'])(7,2,6,:,:,:) = data.eeg.([data.eeg.metrics{n} '_all'])(7,1,6,:,:,:);
          
          % permute to format [roi, roi, freq, more-dims...]:
          data.eeg.(data.eeg.metrics{n}) = permute(data.eeg.([data.eeg.metrics{n} '_all']),[4 5 6 1 2 3]);
        end
        data.eeg.dimName = {'roi','roi','freq','subjEeg','day','cond'};
        data.eeg.dimLabels = {num2cell(1:66),num2cell(1:66),{5,11,25},num2cell(1:10),{1,2},{'pre_co','pre_ce','co','ce','pre_rs','post_rs'}};
        
      else
        if p.eegDatabase==3
          data.eeg = load([paths.databases '/SC_Bastian/icoh_all_lcmvhilbert_20140804.mat']);
          data.eeg.metrics = {'icoh','coh','plv','cohy'};
        else
          if p.eegDatabase==2
            data.eeg = load([paths.databases '/SC_Bastian/icoh_all_mne2_20140730.mat']);
          else
            data.eeg = load([paths.databases '/SC_Bastian/icoh_all_0_20140715.mat']);
          end
          data.eeg.metrics = {'icoh','coh','plv'};
        end
        for n=1:length(data.eeg.metrics)
          % replace NaNs of subject 6 on day 2 with day 1:
          data.eeg.([data.eeg.metrics{n} '_all'])(6,2,:,:,:,:) = data.eeg.([data.eeg.metrics{n} '_all'])(6,1,:,:,:,:);
          
           % permute to format [roi, roi, freq, more-dims...]:
          data.eeg.(data.eeg.metrics{n}) = permute(data.eeg.([data.eeg.metrics{n} '_all']),[5 6 7 1 2 3 4]);
        end
        data.eeg.dimName = {'roi','roi','freq','subjEeg','day','prepost','task'};
        data.eeg.dimLabels = {num2cell(1:66),num2cell(1:66),{5,11,25},num2cell(1:10),{1,2},{'pre','post'},{'co','ce'}};
      end
      data.eeg.dimSize = size(data.eeg.plv);
      
      
      %% pre filter and average over the data:
      modalities=fieldnames(data);
      for k=1:length(modalities)
        myData = data.(modalities{k});
        metrics = myData.metrics;
        
        
        inSpec.dimName = myData.dimName;
        inSpec.dimSize = myData.dimSize;
        inSpec.dimLabels = myData.dimLabels;
        for m=1:length(metrics)
          [myData.(metrics{m}), outSpec] = ConnectomeEnvelopeReduce.filterTensor(myData.(metrics{m}), inSpec, p.(modalities{k}));
        end
        myData.dimName = outSpec.dimName;
        myData.dimSize = outSpec.dimSize;
        myData.dimLabels = outSpec.dimLabels;

%         dimensionNames = fieldnames(p.(modalities{k}));
%         
%         %% calculate indices to filter over all dimensions at once:
%         indices = repmat({':'},size(myData.(metrics{1})));
%         for l=1:length(dimensionNames)
%           curDimId = find(strcmp(myData.dimName,dimensionNames{l}));
%           if isfield(p.(modalities{k}).(dimensionNames{l}),'Ids') && ~isempty(p.(modalities{k}).(dimensionNames{l}).Ids)
%             indices{curDimId} = p.(modalities{k}).(dimensionNames{l}).Ids;
%             myData.dimLabels{curDimId} = myData.dimLabels{curDimId}(indices{curDimId});
%           end
%         end
%         
%         %% now filter all dimensions at once:
%         for m=1:length(metrics)
%           myData.(metrics{m}) = myData.(metrics{m})(indices{:});
%         end
%         
%         %% now compute average over some of the dimensions
%         for l=1:length(dimensionNames)
%           if p.(modalities{k}).(dimensionNames{l}).Avg
%             curDimId = find(strcmp(myData.dimName,dimensionNames{l}));
%             for m=1:length(metrics)
%               myData.(metrics{m}) = mean(myData.(metrics{m}),curDimId);
%             end
%           end
%         end
%         
%         %% squeeze dimensions and at the same time delete corresponding dimNames and dimLabels
%         dimsToReduce = find(size(myData.(metrics{1}))==1);
%         dimsToKeep = find(size(myData.(metrics{1}))~=1);
%         for m=1:length(metrics)
%           myData.(metrics{m}) = permute(myData.(metrics{m}),[dimsToKeep dimsToReduce]);
%         end
%         myData.dimName = myData.dimName(dimsToKeep);
%         myData.dimSize = size(myData.(metrics{1}));
%         myData.dimLabels = myData.dimLabels(dimsToKeep);
%         
        data.(modalities{k}) = myData;
      end
      
      
      %% convert FC-Matrix to Vector:
      numRois = size(data.eeg.plv,1);
      numFreq = size(data.eeg.plv,3);
      trigIds=find(triu(ones(numRois,numRois),1));
      
      for k=1:length(modalities)
        metrics = data.(modalities{k}).metrics;
        for m=1:length(metrics)
          data.(modalities{k}).(metrics{m}) = cell2mat(cellfun(@(x) x(trigIds), num2cell(data.(modalities{k}).(metrics{m}),[1 2]), 'UniformOutput', false));
          data.(modalities{k}).(metrics{m}) = permute(data.(modalities{k}).(metrics{m}),[1 3:length(data.(modalities{k}).dimName) 2]);
        end
        data.(modalities{k}).dimName = cat(2,{'roiPair'},data.(modalities{k}).dimName(3:end));
        data.(modalities{k}).dimSize = size(data.(modalities{k}).(metrics{1}));
        data.(modalities{k}).dimLabels{1} = num2cell(1:length(trigIds));
        data.(modalities{k}).dimLabels(2) = [];
      end
      
      
      
      %% Add Absolute Imaginary Coherence:
      for k=1:length(modalities)
        data.(modalities{k}).aicoh = abs(data.(modalities{k}).icoh);
        data.(modalities{k}).metrics{end+1} = 'aicoh';
      end
      
      %% Start Evaluation:
      metrics = intersect(data.eeg.metrics, data.sim.metrics);
      compareSimExp = struct();
      
      %% 
      if sum(strcmp(data.sim.dimName,'freq'))
        simDimSize = data.sim.dimSize(3:end);
        simDimLabels = data.sim.dimLabels(3:end);
        simDimName = data.sim.dimName(3:end);
      else
        simDimSize = data.sim.dimSize(2:end);
        simDimLabels = data.sim.dimLabels(2:end);
        simDimName = data.sim.dimName(2:end);
      end
      eegDimSize = data.eeg.dimSize(3:end);
      eegDimLabels = data.eeg.dimLabels(3:end);
      eegDimName = data.eeg.dimName(3:end);
      
      
      %% per freq:
      for m=1:length(metrics)
        for freq=1:numFreq
          tmpEEG = reshape(data.eeg.(metrics{m})(:,freq,:,:,:,:,:),[2145 prod(eegDimSize)]);
          if sum(strcmp(data.sim.dimName,'freq'))
            tmpSIM = reshape(data.sim.(metrics{m})(:,freq,:,:,:,:,:),[2145 prod(simDimSize)]);
          else
            tmpSIM = reshape(data.sim.(metrics{m})(:,:,:,:,:,:),[2145 prod(simDimSize)]);
          end
          if strcmp(metrics{m},'cohy')
            compareSimExp.perFreq.(metrics{m}).coh(freq,:,:) = ConnectomeEnvelopeReduce.calcCoherence(tmpEEG',tmpSIM');
          else
            %% Corr:
            [compareSimExp.perFreq.(metrics{m}).rho(freq,:,:), compareSimExp.perFreq.(metrics{m}).pval(freq,:,:)] = corr(tmpEEG,tmpSIM);
            
            %% Jaccard similarity coefficient:
            tmpMin = bsxfun(@min,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
            tmpMax = bsxfun(@max,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
            compareSimExp.perFreq.(metrics{m}).jac(freq,:,:) = sum(tmpMin,1) ./ sum(tmpMax,1);
          end
        end
        if strcmp(metrics{m},'cohy')
          compareSimExp.perFreq.(metrics{m}).coh = reshape(compareSimExp.perFreq.(metrics{m}).coh,[numFreq eegDimSize simDimSize]);
          
          compareSimExp.perFreqAvg.(metrics{m}).coh = permute(mean(compareSimExp.perFreq.(metrics{m}).coh,1),[2:length(size(compareSimExp.perFreq.(metrics{m}).coh)) 1]);
        else
          compareSimExp.perFreq.(metrics{m}).rho = reshape(compareSimExp.perFreq.(metrics{m}).rho,[numFreq eegDimSize simDimSize]);
          compareSimExp.perFreq.(metrics{m}).pval = reshape(compareSimExp.perFreq.(metrics{m}).pval,[numFreq eegDimSize simDimSize]);
          compareSimExp.perFreq.(metrics{m}).jac = reshape(compareSimExp.perFreq.(metrics{m}).jac,[numFreq eegDimSize simDimSize]);
          
          compareSimExp.perFreqAvg.(metrics{m}).rho = permute(mean(compareSimExp.perFreq.(metrics{m}).rho,1),[2:length(size(compareSimExp.perFreq.(metrics{m}).rho)) 1]);
          compareSimExp.perFreqAvg.(metrics{m}).pval = permute(mean(compareSimExp.perFreq.(metrics{m}).pval,1),[2:length(size(compareSimExp.perFreq.(metrics{m}).pval)) 1]);
          compareSimExp.perFreqAvg.(metrics{m}).jac = permute(mean(compareSimExp.perFreq.(metrics{m}).jac,1),[2:length(size(compareSimExp.perFreq.(metrics{m}).jac)) 1]);
        end
      end
      
      %% combine all freq:
      for m=1:length(metrics)
        
        if sum(strcmp(data.sim.dimName,'freq'))
          tmpEEG = reshape(data.eeg.(metrics{m}),[2145*numFreq prod(eegDimSize)]);
          tmpSIM = reshape(data.sim.(metrics{m}),[2145*numFreq prod(simDimSize)]);
        else
          tmpEEG = reshape(mean(data.eeg.(metrics{m}),2),[2145 prod(eegDimSize)]);
          tmpSIM = reshape(data.sim.(metrics{m}),[2145 prod(simDimSize)]);
        end
        
        if strcmp(metrics{m},'cohy')
          compareSimExp.overFreq.(metrics{m}).coh = ConnectomeEnvelopeReduce.calcCoherence(tmpEEG',tmpSIM');
          compareSimExp.overFreq.(metrics{m}).coh = reshape(compareSimExp.overFreq.(metrics{m}).coh,[eegDimSize simDimSize]);
        else
          %% Corr:
          [compareSimExp.overFreq.(metrics{m}).rho, compareSimExp.overFreq.(metrics{m}).pval] = corr(tmpEEG,tmpSIM);
          compareSimExp.overFreq.(metrics{m}).rho = reshape(compareSimExp.overFreq.(metrics{m}).rho,[eegDimSize simDimSize]);
          compareSimExp.overFreq.(metrics{m}).pval = reshape(compareSimExp.overFreq.(metrics{m}).pval,[eegDimSize simDimSize]);
          
          %% Jaccard similarity coefficient:
          tmpMin = bsxfun(@min,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
          tmpMax = bsxfun(@max,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
          compareSimExp.overFreq.(metrics{m}).jac = sum(tmpMin,1) ./ sum(tmpMax,1);
          compareSimExp.overFreq.(metrics{m}).jac = reshape(compareSimExp.overFreq.(metrics{m}).jac, [eegDimSize simDimSize]);
        end
      end
      
      %% combine dimNames and Labels and Sizes:
      compareSimExp.dimSize = [eegDimSize simDimSize];
      compareSimExp.dimLabels = [eegDimLabels simDimLabels];
      compareSimExp.dimName = [eegDimName simDimName];
      
      %% evaluate subject specific:
      eegSubjDim = find(strcmp(compareSimExp.dimName,'subjEeg'));
      simSubjDim = find(strcmp(compareSimExp.dimName,'subjId'));
      if ~isempty(eegSubjDim) && ~isempty(simSubjDim)
        
        nonDiagIds=find(~eye([compareSimExp.dimSize(eegSubjDim) compareSimExp.dimSize(simSubjDim)]));
        
        %% do permutation tests for all metrics:
        dimSize = compareSimExp.dimSize;
        dimSize(end+1:4)=1;
        
        metrics0 = {'overFreq','perFreqAvg'};
        for m0=1:length(metrics0)
          metrics1 = fieldnames(compareSimExp.(metrics0{m0}));
          for m1=1:length(metrics1)
            metrics2 = fieldnames(compareSimExp.(metrics0{m0}).(metrics1{m1}));
            for m2=1:length(metrics2)
              
              if strcmp(metrics2{m2},'rho') || strcmp(metrics2{m2},'jac')
                indices = {':',':',':',':'};
                iterateOverDims = setdiff(1:4,[eegSubjDim simSubjDim]);
                for iterDim1=1:dimSize(iterateOverDims(1))
                  indices{iterateOverDims(1)} = iterDim1;
                  for iterDim2=1:dimSize(iterateOverDims(2))
                    indices{iterateOverDims(2)} = iterDim2;
                    
                    normString = {'orig','norm'};
                    for normBefore = 1:2
                      RhoMatch = cell(2,1);
                      RhoNomat = cell(2,1);
                      
                      transpString = {'fixEeg','fixSim'};
                      for transp=1:2
                        metricVal = atanh(compareSimExp.(metrics0{m0}).(metrics1{m1}).(metrics2{m2})(indices{:}));
                        
                        if transp==2
                          metricVal = metricVal';
                        end
                        
                        if normBefore==2
                          metricVal = bsxfun(@rdivide, metricVal, sum(metricVal,1));
                        end
                        
                        %% extract matching and nonmatching matrix:
                        matching = repmat(diag(metricVal),1,8);
                        ur = triu(metricVal,1);
                        bl = tril(metricVal,-1);
                        if transp==1
                          nonmatch = ur(:,2:end) + bl(:,1:end-1);
                        else
                          nonmatch = (ur(1:end-1,:) + bl(2:end,:))';
                        end
                        
                        %% calc fraction greater vs smaller
                        pairedDiff = matching(:) - nonmatch(:);
                        numMatchSmaller = sum(pairedDiff < 0);
                        numMatchLarger = sum(pairedDiff > 0);
                        fractionMatchLarger = numMatchLarger / (numMatchSmaller+numMatchLarger);
                        compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).fractionMatchLarger(iterDim1,iterDim2) = fractionMatchLarger;
                        
                        %% signrank test
                        if sum(isnan(matching(:))) || sum(isnan(nonmatch(:)))
                          disp(['NaN values in statistical test at m0=' num2str(m0) ' m1=' num2str(m1) ' m2=' num2str(m2) ' iterDim1=' num2str(iterDim1) ' iterDim2=' num2str(iterDim2) ' '  ])
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).signrank_pval(iterDim1,iterDim2) = NaN;
                        else
                          pval = signrank(matching(:),nonmatch(:));
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).signrank_pval(iterDim1,iterDim2) = pval;
                        end
                        
                        %% do paired t-test with repeated nonmatching variables:
                        RhoMatch{transp} = matching(:);
                        RhoNomat{transp} = nonmatch(:);
                        [~,pval]=ttest(RhoMatch{transp},RhoNomat{transp});
                        compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).repeatMatches(iterDim1,iterDim2) = pval;
                        
                        %% do paired t-test where the nonmatching values are before averaged for each column:
                        nonmatchingMeans = mean(nonmatch,2);
                        [~,pval]=ttest(matching(:,1),nonmatchingMeans);
                        compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).meanNonmatch(iterDim1,iterDim2) = pval;
                        
                        %% do paired t-test where the nonmatching values are before averaged for each column:
                        nonmatchingMedian = median(nonmatch,2);
                        [~,pval]=ttest(matching(:,1),nonmatchingMedian);
                        compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).medianNonmatch(iterDim1,iterDim2) = pval;
                        
                        %% mixed-effects anovan:
                        if ~isreal(metricVal(:)) || sum(isnan(matching(:))) || sum(isnan(nonmatch(:)))
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).anovan_pval_matching(iterDim1,iterDim2) = NaN;                        
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).anovan_pval_subjDTI(iterDim1,iterDim2) = NaN;                        
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).anovan_pval_subjEEG(iterDim1,iterDim2) = NaN;                                                  
                        else
                          groupMatch = eye(9,9);                        
                          groupSubjDTI = repmat(1:9,[9 1]);
                          groupSubjEEG = repmat(1:9,[9 1])';
                          pval = anovan(metricVal(:), {groupMatch(:) groupSubjDTI(:) groupSubjEEG(:)} , 'random', [2 3], 'varnames', {'matching', 'subjDTI', 'subjEEG'}, 'display','off');
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).anovan_pval_matching(iterDim1,iterDim2) = pval(1);                        
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).anovan_pval_subjDTI(iterDim1,iterDim2) = pval(2);                        
                          compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).(transpString{transp}).(normString{normBefore}).anovan_pval_subjEEG(iterDim1,iterDim2) = pval(3);                        
                        end
                        
                      end
                      
                      [~,compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).fixBoth(iterDim1,iterDim2)]=ttest([RhoMatch{1}(:); RhoMatch{2}(:)],[RhoNomat{1}(:); RhoNomat{2}(:)]);
                      
                    end
                    
                  end
                end

                iterateOverDims(iterateOverDims>length(compareSimExp.dimName)) = [];
                compareSimExp.permTest.dimName = compareSimExp.dimName(iterateOverDims);
                compareSimExp.permTest.dimSize = size(compareSimExp.permTest.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}).fixSim);
                compareSimExp.permTest.dimLabels = compareSimExp.dimLabels(iterateOverDims);
              end
            end
          end
        end
      end
      
      %% maybe combine subject dimensions from eeg and sim into one dimension:
      if ~isempty(eegSubjDim) && ~isempty(simSubjDim)
        if strcmp(p.resultsCombineSubjDims,'match') || strcmp(p.resultsCombineSubjDims,'nonmatch')
          nonDiagIds=find(~eye([compareSimExp.dimSize(eegSubjDim) compareSimExp.dimSize(simSubjDim)]));
          
          metrics0 = {'overFreq','perFreqAvg'};
          for m0=1:length(metrics0)
            metrics1 = fieldnames(compareSimExp.(metrics0{m0}));
            for m1=1:length(metrics1)
              metrics2 = fieldnames(compareSimExp.(metrics0{m0}).(metrics1{m1}));
              for m2=1:length(metrics2)                
                asCell = num2cell(compareSimExp.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}),[eegSubjDim simSubjDim]);
                asCell = cellfun(@(x) squeeze(x), asCell,'UniformOutput', false);
                if strcmp(p.resultsCombineSubjDims,'match')
                  asCell = cellfun(@(x) shiftdim(diag(x),eegSubjDim-1), asCell,'UniformOutput', false);
                elseif strcmp(p.resultsCombineSubjDims,'nonmatch')
                  asCell = cellfun(@(x) shiftdim(x(nonDiagIds),eegSubjDim-1), asCell,'UniformOutput', false); %#ok<FNDSB>
                end
                compareSimExp.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}) = permute(cell2mat(asCell),[setdiff(1:length(compareSimExp.dimSize),simSubjDim) simSubjDim]);
                
                % calc best param per subject
                tmp = compareSimExp.(metrics0{m0}).(metrics1{m1}).(metrics2{m2});
                siz = size(tmp);
                numParam = length(siz);
                siz = siz(2:end);
                siz(end+1:5) = 1;
                [~,bestLinInd] = max(reshape(tmp,[9 prod(siz)]),[],2);
                [I1,I2,I3,I4,I5] = ind2sub(siz,bestLinInd);
                allPerSubj = [I1,I2,I3,I4,I5];
                compareSimExp.bestParamPerSubj.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}) = allPerSubj(:,1:numParam-1);
                
              end
            end
          end
          
          compareSimExp.dimSize = size(compareSimExp.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}));
          compareSimExp.dimLabels{eegSubjDim} = num2cell(1:compareSimExp.dimSize(eegSubjDim));
          compareSimExp.dimName(simSubjDim) = [];
          compareSimExp.dimLabels(simSubjDim) = [];
        end
        
      end
      
      %% maybe reduce some dimensions in the results:
      if isfield(p,'results')
        inSpec.dimName = compareSimExp.dimName;
        inSpec.dimSize = compareSimExp.dimSize;
        inSpec.dimLabels = compareSimExp.dimLabels;
        metrics0 = {'overFreq','perFreqAvg'};
        for m0=1:length(metrics0)
          metrics1 = fieldnames(compareSimExp.(metrics0{m0}));
          for m1=1:length(metrics1)
            metrics2 = fieldnames(compareSimExp.(metrics0{m0}).(metrics1{m1}));
            for m2=1:length(metrics2)

                [compareSimExp.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}), outSpec] = ...
                  ConnectomeEnvelopeReduce.filterTensor(...
                  compareSimExp.(metrics0{m0}).(metrics1{m1}).(metrics2{m2}), ...
                  inSpec, ...
                  p.results);

            end
          end
        end
        compareSimExp.dimName = outSpec.dimName;
        compareSimExp.dimSize = outSpec.dimSize;
        compareSimExp.dimLabels = outSpec.dimLabels;
      end
      
      %% calc the following only for the simulation with the best coh value:
      if sum(strcmp(metrics,'cohy')) && sum(strcmp(data.sim.dimName,'freq'))
        [~,Ibest] = max(compareSimExp.overFreq.coh.rho(:));
        [compareSimExp.samples.bestEegId,compareSimExp.samples.bestSimId] = ind2sub([prod(eegDimSize) prod(simDimSize)],Ibest);
        tmpEEG = reshape(data.eeg.cohy,[2145*numFreq prod(eegDimSize)]);
        tmpSIM = reshape(data.sim.cohy,[2145*numFreq prod(simDimSize)]);
        compareSimExp.samples.cohy_eegBest = tmpEEG(:,compareSimExp.samples.bestEegId);
        compareSimExp.samples.cohy_simBest = tmpSIM(:,compareSimExp.samples.bestSimId);
      end
      
      compareSimExp.dimName_eeg = data.eeg.dimName;
      compareSimExp.dimName_sim = data.sim.dimName;
      compareSimExp.dims_eeg = data.eeg.dimSize;
      compareSimExp.dims_sim = data.sim.dimSize;
      compareSimExp.dimLabels_eeg = data.eeg.dimLabels;
      compareSimExp.dimLabels_sim = data.sim.dimLabels;
      
      compareSimExp.sim.dimName = data.sim.paramName;
      compareSimExp.sim.dimLabels = data.sim.paramVar;
      if isfield(data.sim,'meanOrderParam')
        compareSimExp.sim.metrics.orderParamMean = data.sim.meanOrderParam;
        compareSimExp.sim.metrics.orderParamStd = data.sim.stdOrderParam;
      end
      
      save(fullfile(savepath,['compareSimExp' num2str(this.currJobid)]),'-struct','compareSimExp')
    end
    
    
    
    function plotResults(this,compareSimExp)
      
      p = this.params.ConnectomeEnvelopeReduce;
      savepath = fullfile(this.workpath,p.outDirectory);
      plotDir = fullfile(savepath,['plots' num2str(this.currJobid)]);
      mkdir(plotDir)
      
      if nargin<2
        compareSimExp = load(fullfile(savepath,['compareSimExp' num2str(this.currJobid)]));
      end
      
      dimNames = compareSimExp.dimName;
      dimSizes = compareSimExp.dimSize;
      dimLabels = compareSimExp.dimLabels;
      
      
      %% plot normal result metrics:
      ConnectomeEnvelopeReduce.plotStructure(compareSimExp.overFreq, compareSimExp, 'result.overFreq', 'result_overFreq', plotDir, {'pval'}, p.permutePlotDims)
      ConnectomeEnvelopeReduce.plotStructure(compareSimExp.perFreqAvg, compareSimExp, 'result.perFreqAvg', 'result_perFreqAvg', plotDir, {'pval'}, p.permutePlotDims)
      
      %% plot over fixed subject if two remaining dimensions correspond to eeg and mri subject
      if length(compareSimExp.dimName)==2 && sum(strcmp(compareSimExp.dimName,'subjEeg')) && sum(strcmp(compareSimExp.dimName,'subjId'))
        ConnectomeEnvelopeReduce.plotStructureSubj(compareSimExp.overFreq, compareSimExp, 'fixSubj.overFreq', 'fixSubj_overFreq', plotDir, {'pval'})
        ConnectomeEnvelopeReduce.plotStructureSubj(compareSimExp.perFreqAvg, compareSimExp, 'fixSubj.perFreqAvg', 'fixSubj_perFreqAvg', plotDir, {'pval'})
      end
      
      
      %% plot permutation eval:
      if isfield(compareSimExp,'permTest') && p.plotPermTests
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.permTest.overFreq, compareSimExp.permTest, 'ttest.overFreq', 'ttest_overFreq', plotDir, {});
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.permTest.perFreqAvg, compareSimExp.permTest, 'ttest.perFreqAvg', 'ttest_perFreqAvg', plotDir, {});
      end
      
      %% write permutation eval to file:
      if isfield(compareSimExp,'permTest') && isempty(compareSimExp.permTest.dimName)
        fileID = fopen(fullfile(plotDir,'subjPermTest_pval.txt'),'w');
        ConnectomeEnvelopeReduce.dispStructure(compareSimExp.permTest.overFreq, 'permTest.overFreq', fileID, {});
        ConnectomeEnvelopeReduce.dispStructure(compareSimExp.permTest.perFreqAvg, 'permTest.perFreqAvg', fileID, {});
        fclose(fileID);
      end
      
      
      %% plot coh of cohy
      if isfield(compareSimExp,'samples') && isfield(compareSimExp.samples,'cohy_eegBest') && isfield(compareSimExp.samples,'cohy_simBest')
        cohy_eegBest = compareSimExp.samples.cohy_eegBest;
        cohy_simBest = compareSimExp.samples.cohy_simBest;
        
        %% calc desc:
        if length(compareSimExp.dims_eeg)==4
          [I(1),I(2)] = ind2sub(compareSimExp.dims_eeg(3:end),compareSimExp.samples.bestEegId);
        elseif length(compareSimExp.dims_sim)==3
          I(1) = compareSimExp.samples.bestEegId;
          I(2) = compareSimExp.samples.bestSimId;
        else
          [I(1),I(2),I(3),I(4)] = ind2sub(compareSimExp.dims_sim(3:end),compareSimExp.samples.bestSimId);
        end
        descParam = '';
        for ii=1:min(length(dimNames),length(I))
          descParam = [descParam dimNames{ii} '=' num2str(dimLabels{ii}{I(ii)}) ' ']; %#ok<AGROW>
        end
        
        %% plots:
        figure(1); clf;
        bothCoh=abs(cohy_simBest).*abs(cohy_eegBest);
        phaseDiffDirect=mod(angle(cohy_simBest)-angle(cohy_eegBest)+pi,2*pi)-pi;
        scatter(bothCoh,phaseDiffDirect)
        title(['Phasediff vs coh for ' descParam])
        xlabel('eegCoh * simCoh')
        ylabel('eeg-sim-difference in phase-lag between rois [rad]')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        print('-dpng','-r72',fullfile(plotDir,'cohy_phasediff_bestSim.png'))
        
        figure(1); clf;
        bothCoh=abs(cohy_simBest).*abs(cohy_eegBest);
        phaseDiffDirect=mod(angle(cohy_simBest)-angle(cohy_eegBest)+pi,2*pi)-pi;
        scatter(bothCoh,abs(phaseDiffDirect))
        title(['Phasediff vs coh for ' descParam])
        xlabel('eegCoh * simCoh')
        ylabel('abs eeg-sim-difference in phase-lag between rois [rad]')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        print('-dpng','-r72',fullfile(plotDir,'cohy_absPhasediff_bestSim.png'))
        
        figure(1); clf;
        hist(angle(cohy_eegBest),100)
        title(['Histogram over all ROI-pairs eeg'])
        xlabel('phase lag [rad]')
        ylabel('number of ROI-pairs')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        print('-dpng','-r72',fullfile(plotDir,'hist_phaseLag_EEG.png'))
        
        figure(1); clf;
        hist(angle(cohy_simBest),100)
        title(['Histogram over all ROI-pairs for sim'])
        xlabel('phase lag [rad]')
        ylabel('number of ROI-pairs')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        print('-dpng','-r72',fullfile(plotDir,'hist_phaseLag_bestSim.png'))
      end
      
      %% plot sim specifics
      if isfield(compareSimExp,'sim')
        if isfield(compareSimExp.sim,'metrics')
          metrics = fieldnames(compareSimExp.sim.metrics);
          for m=1:length(metrics)
            figure(1); clf;
            imagesc(compareSimExp.sim.metrics.(metrics{m})(:,:,1,1));
            title(metrics{m})
            colormap hot;
            colorbar
            if length(compareSimExp.sim.dimName)>1
              xlabel(compareSimExp.sim.dimName{2})
              set(gca,'XTick',1:length(compareSimExp.sim.dimLabels{2}));
              set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%10.1f'),compareSimExp.sim.dimLabels{2},'UniformOutput',false))
            end
            ylabel(compareSimExp.sim.dimName{1})
            set(gca,'YTick',1:length(compareSimExp.sim.dimLabels{1}));
            set(gca,'YTickLabel',compareSimExp.sim.dimLabels{1})
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
            print('-dpng','-r72',fullfile(plotDir,[metrics{m} '.png']))
          end
        end
      end
      
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
  
  methods(Static)
    function coh = calcCoherence(sig1,sig2)
      autospectrum_eeg = mean( conj(sig1) .* sig1 , 2 );
      autospectrum_sim = mean( conj(sig2) .* sig2 , 2 );
      crossspectrum = permute(mean( bsxfun(@times, conj(sig1), permute(sig2,[3 2 1])) , 2 ),[1 3 2]);
      coh = abs(crossspectrum ./ sqrt( autospectrum_eeg * autospectrum_sim' ));
    end
    
    function dispStructure(structure, caption, fileID, excludeFields)
      % use as follows:
      % fileID = fopen('permTestResults.txt','w');
      % ConnectomeEnvelopeReduce.dispStructure(structure, 'permTest', fileID, {})
      % fclose(fileID);
      
      if isstruct(structure)
        %% recursive call:
        fnames = fieldnames(structure);
        for f=1:length(fnames)
          if sum(strcmp(excludeFields,fnames{f}))==0
            ConnectomeEnvelopeReduce.dispStructure(structure.(fnames{f}), [caption '.' fnames{f}], fileID, excludeFields);
          end
        end
      else
        fprintf(fileID,'%s = %12.8f\n',caption,structure);
      end
      
    end
    
    
    function plotStructure(structure, inSpec, caption, filename, plotDir, excludeFields, permutePlotDims)
      
      if nargin<7
        permutePlotDims = [];
      end
      
      if isstruct(structure)
        
        %% recursive call:
        fnames = fieldnames(structure);
        for f=1:length(fnames)
          if sum(strcmp(excludeFields,fnames{f}))==0
            ConnectomeEnvelopeReduce.plotStructure(structure.(fnames{f}), inSpec, [caption '.' fnames{f}], [filename '_' fnames{f}], plotDir, excludeFields, permutePlotDims);
          end
        end
        
      else
        
        %% plot:
        siz = size(structure);
        
        if ~isempty(permutePlotDims)
          inSpecPerm = inSpec;
          structure = permute(structure,permutePlotDims);
          while length(inSpecPerm.dimName)<length(permutePlotDims)
            inSpecPerm.dimName{end+1} = '';
            inSpecPerm.dimLabels{end+1} = {};
          end
          inSpecPerm.dimName = inSpecPerm.dimName(permutePlotDims);
          inSpecPerm.dimLabels = inSpecPerm.dimLabels(permutePlotDims);
        else
          inSpecPerm = inSpec;
        end
        
        if length(siz) > 2
          
%           if length(siz) == 3
%             
%             m=round(sqrt(siz(3)));
%             n=ceil(siz(3)/m);
%             
%             figure(1); clf;
%             for k=1:siz(3)
%               subplot(m,n,k);
%               imagesc(structure(:,:,k));
%               title(caption)
%               colormap hot;
%               colorbar;
%               if length(inSpecPerm.dimName)>1
%                 xlabel(inSpecPerm.dimName{2})
%                 set(gca,'XTick',1:size(structure,2));
%                 set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),inSpecPerm.dimLabels{2},'UniformOutput',false))
%               end
%               ylabel(inSpecPerm.dimName{1})
%               set(gca,'YTick',1:size(structure,1));
%               set(gca,'YTickLabel',inSpecPerm.dimLabels{1})
%             end
%             set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
%             print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
%             
%           elseif length(siz) == 4
            
            minc = min(structure(:));
            maxc = max(structure(:));
            
            figure(1)
            clf
            p = panel();
            p.pack(size(structure,3), size(structure,4));
            p.de.margin = 2;
            p.margin = [25 25 20 2];
            p.fontsize = 8;
            for m = 1:size(structure,3)
              for n = 1:size(structure,4)
                p(m, n).select();
                imagesc(structure(:,:,m,n));
                colormap hot;
                set(gca,'clim',[minc maxc]);
                axis image;
                
                if n==size(structure,4)
                  originalSize = get(gca, 'Position');
                  colorbar;
                  set(gca, 'Position', originalSize);
                end
                
                if m==size(structure,3)
                  if n==round(size(structure,4)/2)
                    %% with outer xlabel
                    p(m,n).xlabel({inSpecPerm.dimName{2},[],['{\bf' num2str(inSpecPerm.dimLabels{4}{n}) '}'],['{\bf' inSpecPerm.dimName{4} '}']});
                  else
                    %% only inner xlabel
                    p(m,n).xlabel({inSpecPerm.dimName{2},[],['{\bf' num2str(inSpecPerm.dimLabels{4}{n}) '}']});
                  end
                  %set(gca,'XTick',1:size(structure,1));
                  xticks = get(gca,'XTick');
                  set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),inSpecPerm.dimLabels{2}(xticks),'UniformOutput',false))
                  %set(gca,'XTickMode','auto')
                else
                  set(gca, 'xticklabel', {});
                end
                if n==1
                  if isempty(inSpecPerm.dimLabels{3})
                    p(m,n).ylabel(inSpecPerm.dimName{1});
                  else
                    if m==round(size(structure,3)/2)
                      %% with outer ylabel
                      p(m,n).ylabel({['{\bf' inSpecPerm.dimName{3} '}'],['{\bf' num2str(inSpecPerm.dimLabels{3}{m}) '}'],[],inSpecPerm.dimName{1}});
                    else
                      %% only inner ylabel
                      p(m,n).ylabel({['{\bf' num2str(inSpecPerm.dimLabels{3}{m}) '}'],[],inSpecPerm.dimName{1}});
                    end
                  end
                  yticks = get(gca,'YTick');
                  set(gca,'YTickLabel',inSpecPerm.dimLabels{1}(yticks))
                else
                  set(gca, 'yticklabel', {});
                end
                
              end
            end
            
%             set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
%             print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
%             export_fig(fullfile(plotDir,[filename '.pdf']))
%             p.export(fullfile(plotDir,[filename '.pdf']), '-a1.4', '-rp');
            p.export(fullfile(plotDir,[filename '.pdf']), '-rp','-w300','-h100');
            
%           else
%             return;
%           end
        elseif siz(2) == 1 && siz(1) == 1
          %% return because only scalar
          return;
        elseif siz(2) == 1
          figure(1); clf;
          plot(structure);
          title(caption)
          xlabel(inSpecPerm.dimName{1})
          ylabel(caption)
          set(gca,'XTick',1:size(structure,1));
          set(gca,'XTickLabel',inSpecPerm.dimLabels{1})
          set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
          print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
        else
          figure(1);clf;
          imagesc(structure);
          title(caption)
          colormap hot;
          colorbar;
          if length(inSpecPerm.dimName)>1
            xlabel(inSpecPerm.dimName{2})
            set(gca,'XTick',1:size(structure,2));
            set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),inSpecPerm.dimLabels{2},'UniformOutput',false))
          end
          ylabel(inSpecPerm.dimName{1})
          set(gca,'YTick',1:size(structure,1));
          set(gca,'YTickLabel',inSpecPerm.dimLabels{1})
          set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
          print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
        end
        
      end
      
    end
    
    function plotStructureSubj(structure, inSpec, caption, filename, plotDir, excludeFields)
      
      if isstruct(structure)
        
        %% recursive call:
        fnames = fieldnames(structure);
        for f=1:length(fnames)
          if sum(strcmp(excludeFields,fnames{f}))==0
            ConnectomeEnvelopeReduce.plotStructureSubj(structure.(fnames{f}), inSpec, [caption '.' fnames{f}], [filename '_' fnames{f}], plotDir, excludeFields);
          end
        end
        
      else

        %% plot eeg-subj vs mri-subj:
        ur = triu(structure,1);
        bl = tril(structure,-1);
        match = diag(structure);
        fixSubj = {'eeg','mri'};
        for f=1:length(fixSubj)
          if strcmp(fixSubj{f},'eeg')
            nonmatch = ur(:,2:end) + bl(:,1:end-1);
          else
            nonmatch = (ur(1:end-1,:) + bl(2:end,:))';
          end

          %% plot matching and nonmatching data
          figure(1); clf
          plot(match,'bo')
          hold on;
          plot(mean(nonmatch,2),'g*')
          plot(median(nonmatch,2),'gx')
          plot(nonmatch,'r.','MarkerSize', 12)
          hold off;
          xlabel([fixSubj{f} ' subj'])
          set(gca,'XTick',1:length(match));
          if strcmp(fixSubj{f},'eeg')
            set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),inSpec.dimLabels{1},'UniformOutput',false))
          else
            set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),inSpec.dimLabels{2},'UniformOutput',false))
          end
          ylabel(caption)
          title(caption)
          xlim([0 10])
          legend('match','mean(nonmatch)','median(nonmatch)','nonmatch')
          set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
          print('-dpng','-r72',fullfile(plotDir,[filename '_fix' fixSubj{f} '.png']))
          
        end
      
      end
      
    end
      
    
    
    function [tensor, outSpec] = filterTensor(tensor, inSpec, filterSpec)
%       tensor = rand([5,10,2]);
%       filterSpec.subj.Ids = [1:4 7:10];
%       filterSpec.subj.Avg = false;
%       filterSpec.day.Ids = [];
%       filterSpec.day.Avg = true;
%       inSpec.dimName = {'day', 'subjEeg', 'cond'};
%       inSpec.dimSize = [5 10 2];
%       inSpec.dimLabels = {'', '', ''};
      
      filtDimNames = fieldnames(filterSpec);
      
      %% calculate indices to filter over all dimensions at once:
      indices = repmat({':'},size(inSpec.dimSize));
      for l=1:length(filtDimNames)
        curDimId = find(strcmp(inSpec.dimName,filtDimNames{l}));
        if isfield(filterSpec.(filtDimNames{l}),'Ids') && ~isempty(filterSpec.(filtDimNames{l}).Ids)
          indices{curDimId} = filterSpec.(filtDimNames{l}).Ids;
          inSpec.dimLabels{curDimId} = inSpec.dimLabels{curDimId}(indices{curDimId});
        end
      end
      
      %% now filter all dimensions at once:
      tensor = tensor(indices{:});
      
      %% now compute average over some of the dimensions
      for l=1:length(filtDimNames)
        if filterSpec.(filtDimNames{l}).Avg
          curDimId = find(strcmp(inSpec.dimName,filtDimNames{l}));
          tensor = mean(tensor,curDimId);
        end
      end
      
      %% squeeze dimensions and at the same time delete corresponding dimNames and dimLabels
      dimsToReduce = find(size(tensor)==1);
      dimsToKeep = find(size(tensor)~=1);
      tensor = permute(tensor,[dimsToKeep dimsToReduce]);
      outSpec.dimName = inSpec.dimName(dimsToKeep);
      outSpec.dimSize = size(tensor);
      outSpec.dimLabels = inSpec.dimLabels(dimsToKeep);
      
    end
    
  end
  
end

