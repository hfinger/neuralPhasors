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
      this.params.ConnectomeEnvelopeReduce.compareSC = false;
      this.params.ConnectomeEnvelopeReduce.compareSC_shuffle = false;
      this.params.ConnectomeEnvelopeReduce.compareSC_shufflePermutations = [];
      this.params.ConnectomeEnvelopeReduce.compareSC_db = 1;
      this.params.ConnectomeEnvelopeReduce.useEnvFreqAsParamVar = false;
      
      this.params.ConnectomeEnvelopeReduce.eegDatabase = 1; 
      % 1=icoh_all_0_20140715.mat initial lcmv
      % 2=icoh_all_mne2_20140730.mat uses mne
      % 3=icoh_all_lcmvhilbert_20140804.mat includes cohy uses lcmv and hilbert
      % 4=icoh_all_lcmvhilbertrest_20140807.mat only real RS (also lcmv)
      % 5=icoh_all_lcmv0_001hilbert_rest_bp_20140902.mat
      % 6=icoh_all_lcmv0_001_bp_hilbert_20140907.mat
      % 7=icoh_eeg_funcconn_bp_lcmv0_001_hilbert_20140912.mat
      % 8=icoh_patients_eeg_funcconn_bp_lcmv0001_hilbert_20140912.mat
      % 9=patients_eeg_funcconn_bp_lcmv0001_hilbert_20140930.mat
      % 10=icoh_eeg_funcconn_bp_lcmv0_001_hilbert_20141002.mat (replaces 7)
      % 11=eeg_20141002_funcconn_bp_elor0_001_hilbert.mat
      % 12=eeg_20141003_funcconn_bp_mne_prewhit_0_001_hilbert.mat
      % 13=eeg_20141003_funcconn_bp_mne_noprewhit_0_001_hilbert.mat
      % 14=eeg_20141008_funcconn_bp_lcmv0_001_hilbert_LPC.mat
      % 15=20141006_conn_bp.mat this is in electrode space
      % 16=eeg_20141028_funconn_bp11_lcmv_hilbert.mat includes wPLI and PLI
      % 17=eeg_20141103_controls_fs_funconn_bp_lcmv_hilbert_3_30.mat includes wPLI and PLI
      % 18=eeg_20141109_controls_fs_funconn_bp_pinvonall_hilbert_3_30.mat
      % 19=eeg_20141109_controls_fs_funconn_bp_pinvperroi_hilbert_3_30.mat
      % 20=eeg_20141212_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30.mat
      % 21=eeg_20150102_controls_fs_funconn_bponetrial_elor_hilbert_51122_020115.mat
      % 22=eeg_20141124_controls_bponetrial_elor_hilbert_51122.mat
      % 23=eeg_20150113_controls_fs_funconn_lcmv_bponetrial_hilbert_3_30.mat
      
      % 24=eeg_20150114_controls_fs_funconn_elor_bponetrial_hilbert_51122_120115.mat % from here fixed mirroring of FC
      % 25=eeg_20150114_controls_fs_funconn_lcmv_bponetrial_hilbert_3_30.mat
      % 26=eeg_20150114_conn_bponetrial_3_30.mat
      % 27=eeg_20150125_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30_entmirrored.mat
      % 28=eeg_20150206_fc_patients_mne_bp_hilbert.mat
      % 29=eeg_20150208_fc_patients_mne0_001_bp_hilbert.mat
      
      
      
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
      this.params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder = []; % if empty use jobid
      this.params.ConnectomeEnvelopeReduce.savePreEvaluationData = false;
      
      this.params.ConnectomeEnvelopeReduce.resultsCombineSubjDims = 'no'; % 'no' or 'match' or 'nonmatch' or 'avgMatch'
      
      this.params.ConnectomeEnvelopeReduce.results = struct();
%       this.params.ConnectomeEnvelopeReduce.results.Dim1.Ids = [];
%       this.params.ConnectomeEnvelopeReduce.results.Dim1.Avg = false;
      this.params.ConnectomeEnvelopeReduce.calcSubjectSpecificTests = false;
      this.params.ConnectomeEnvelopeReduce.calcSquaredDist = false;
      this.params.ConnectomeEnvelopeReduce.calcSquaredDistAvg = false;
      this.params.ConnectomeEnvelopeReduce.calcPartialCorrEuclDist = false;
      this.params.ConnectomeEnvelopeReduce.calcPartialCorrFiberDist = false;
      this.params.ConnectomeEnvelopeReduce.calcCrossFreq = false;
      this.params.ConnectomeEnvelopeReduce.calcPhaseLags = false;
      this.params.ConnectomeEnvelopeReduce.calcCompareAicohWithCoh = false;
      this.params.ConnectomeEnvelopeReduce.calcFisher = false;
      this.params.ConnectomeEnvelopeReduce.doNormalizeSubjectMat = false;
      
      this.params.ConnectomeEnvelopeReduce.doPlot = true; 
      this.params.ConnectomeEnvelopeReduce.deleteEnvelopeWhenDone = false;
      this.params.ConnectomeEnvelopeReduce.deleteSimWhenDone = false;
      this.params.ConnectomeEnvelopeReduce.permutePlotDims = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsPerFreq = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsPermTest = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsSquaredDist = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsSquaredDistAvg = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsCrossFreq = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsPhaseLags2DHist = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsPhaseLagsEval = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsPhaseBinsEeg = [];
      this.params.ConnectomeEnvelopeReduce.permutePlotDimsPhaseBinsSim = [];
      this.params.ConnectomeEnvelopeReduce.plotPermTests = true;
      this.params.ConnectomeEnvelopeReduce.plotBestParamPerSubj = false;
      this.params.ConnectomeEnvelopeReduce.excludePlotFieldsRegexp = '';
      this.params.ConnectomeEnvelopeReduce.deletePlotFolder = false;
      
      this.params.ConnectomeEnvelopeReduce.onlyCollectConnFC = false;
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
      if p.onlyCollectConnFC
        
        disp('Start calculating ConnFC...')
        this.gatherResults();
        
      else
      
        if p.reloadCompareSimExp
          disp('Start loading compareSimExp...')
          compareSimExp = load(fullfile(savepath,['compareSimExp' num2str(this.currJobid)]));
        else

          if ischar(p.reloadConnFC)
            ConnFC = load(fullfile(savepath,p.reloadConnFC));
          elseif p.reloadConnFC
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
      
      if p.compareSC
        paths = dataPaths( );
        if p.compareSC_db==1 || p.compareSC_db==4
          
          if p.compareSC_db==4
            dataSCTmp = load([paths.databases '/SC_Bastian/dti_20141209_preprocessed.mat']);
          elseif p.compareSC_db==1
            dataSCTmp = load([paths.databases '/SC_Bastian/dist_and_CI_controls_preprocessed.mat']);
          end
          
          numFreq = 1;
          numRois = 66;
          
          if p.compareSC_shuffle
            trigIds = find(triu(ones(size(dataSCTmp.ci{1})),1));
            
            if ~isempty(p.compareSC_shufflePermutations)
              permuteTmp = p.compareSC_shufflePermutations;
            else
              permuteTmp = randperm(length(trigIds(:)),length(trigIds(:)));
            end
            
            for subj=1:length(dataSCTmp.ci)
              if ~isempty(dataSCTmp.ci{subj})
                tmp = dataSCTmp.ci{subj}(trigIds);
                tmp = tmp(permuteTmp);
                dataSCTmp.ci{subj} = zeros(size(dataSCTmp.ci{1}));
                dataSCTmp.ci{subj}(trigIds) = tmp;
              end
            end
          end
          
          SCs = cell2mat(permute(dataSCTmp.ci,[1 3 4 2]));
          SCs(isnan(SCs)) = 0;
          SCs = SCs + permute(SCs,[2 1 3 4]);
          
          plv = SCs;
          icoh = SCs;
          coh = SCs;
          cohy = SCs;
          FCsim = SCs;
          lpc = SCs;
          pli = SCs;
          wpli = SCs;
          
          paramName{1} = 'subjId';
          if p.compareSC_db==4
            paramVar{1} = num2cell([1:4 6:13 15 17:22]);
          elseif p.compareSC_db==1
            paramVar{1} = num2cell([1:4 6:10]);
          end
          dims = length(paramVar{1});
        else
          paths = dataPaths( );
          dataSCTmp = load([paths.databases '/SC_Bastian/patients_t1_logCI_mul_20140924_preprocessed.mat']);
          
          numFreq = 1;
          numRois = 66;
          
          SCs = cell2mat(permute(dataSCTmp.ci,[1 3 4 2]));
          SCs(isnan(SCs)) = 0;
          SCs = SCs + permute(SCs,[2 1 3 4]);
          
          plv = SCs;
          icoh = SCs;
          coh = SCs;
          cohy = SCs;
          FCsim = SCs;
          lpc = SCs;
          pli = SCs;
          wpli = SCs;
          
          paramName{1} = 'subjId';
          paramVar{1} = num2cell([7    12    14    15    16    22    24    25]);
          dims = 8;
        end
        
      else
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
            lpc = zeros(numRois,numRois,numFreq,numSims);
            pli = zeros(numRois,numRois,numFreq,numSims);
            wpli = zeros(numRois,numRois,numFreq,numSims);
            for k=1:numSims
              stats=load([this.workpath '/' p.ConnectomeSimOut '/' num2str(k) 'FC.mat']);

              plv(:,:,1,k) = stats.FCsimNoBold;
              icoh(:,:,1,k) = stats.FCsimNoBold;
              coh(:,:,1,k) = stats.FCsimNoBold;
              cohy(:,:,1,k) = stats.FCsimNoBold;
              FCsim(:,:,1,k) = stats.FCsimNoBold;
              lpc(:,:,1,k) = stats.FCsimNoBold;
              pli(:,:,1,k) = stats.FCsimNoBold;
              wpli(:,:,1,k) = stats.FCsimNoBold;
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
          lpc = zeros(numRois,numRois,numFreq,numSims);
          pli = zeros(numRois,numRois,numFreq,numSims);
          wpli = zeros(numRois,numRois,numFreq,numSims);
          meanOrderParam = zeros(1,numSims);
          stdOrderParam = zeros(1,numSims);
          for k=1:numSims
%             disp(k)
            stats=load([this.workpath '/' p.ConnectomeEnvelopeOut '/FCsimJob' num2str(k) '.mat']);
            for freq=1:numFreq
              plv(:,:,freq,k) = stats.PLVsim{freq}.PLV{1};
              icoh(:,:,freq,k) = stats.ICOHsim{freq}.icoh{1};
              coh(:,:,freq,k) = stats.ICOHsim{freq}.coherence{1};
              cohy(:,:,freq,k) = stats.ICOHsim{freq}.coherency{1};
              FCsim(:,:,freq,k) = stats.FCsim{freq}.FC{1};
              lpc(:,:,freq,k) = imag(stats.ICOHsim{freq}.coherency{1})./sqrt(1-real(stats.ICOHsim{freq}.coherency{1}).^2);
              pli(:,:,freq,k) = stats.ICOHsim{freq}.pli{1};
              wpli(:,:,freq,k) = stats.ICOHsim{freq}.wpli{1};
            end

            if exist(fullfile(this.workpath, p.ConnectomeSimOut),'dir')
              statsFile = [this.workpath '/' p.ConnectomeSimOut '/' num2str(k) 'KuramotoStats.mat'];
              if exist(statsFile,'file')
                simStats=load(statsFile);
                meanOrderParam(k) = simStats.meanOrderParam;
                stdOrderParam(k) = simStats.stdOrderParam;
              end
            end
          end


          if p.useEnvFreqAsParamVar

            plv = permute(plv,[1 2 5 4 3]);
            icoh = permute(icoh,[1 2 5 4 3]);
            coh = permute(coh,[1 2 5 4 3]);
            cohy = permute(cohy,[1 2 5 4 3]);
            FCsim = permute(FCsim,[1 2 5 4 3]);
            lpc = permute(lpc,[1 2 5 4 3]);
            pli = permute(pli,[1 2 5 4 3]);
            wpli = permute(wpli,[1 2 5 4 3]);

            dims(end+1) = numFreq;
            paramName{end+1} = 'envFreq';
            paramVar{end+1} = num2cell(cellfun(@(x) (x.bp.Fp1+x.bp.Fp2)/2, stats.ICOHsim));

            numFreq = 1;
          end


        end
      
      end
      plv = reshape(plv,[numRois numRois numFreq dims]);
      icoh = reshape(icoh,[numRois numRois numFreq dims]);
      coh = reshape(coh,[numRois numRois numFreq dims]);
      cohy = reshape(cohy,[numRois numRois numFreq dims]);
      FCsim = reshape(FCsim,[numRois numRois numFreq dims]);
      lpc = reshape(lpc,[numRois numRois numFreq dims]);
      pli = reshape(pli,[numRois numRois numFreq dims]);
      wpli = reshape(wpli,[numRois numRois numFreq dims]);
      
      if length(dims)==1
        dims=[dims 1];
      end
      
      if exist('meanOrderParam','var') && ~p.useEnvFreqAsParamVar
        meanOrderParam = reshape(meanOrderParam,dims);
        stdOrderParam = reshape(stdOrderParam,dims);
        ConnFC.meanOrderParam = meanOrderParam;
        ConnFC.stdOrderParam = stdOrderParam;
      end
      
      ConnFC.sim.plv = plv;
      ConnFC.sim.icoh = icoh;
      ConnFC.sim.coh = coh;
      ConnFC.sim.cohy = cohy;
      ConnFC.sim.FCsim = FCsim;
      ConnFC.sim.lpc = lpc;
      ConnFC.sim.pli = pli;
      ConnFC.sim.wpli = wpli;
      
      ConnFC.sim.spec.dimName = cat(2,{'roi','roi','freq'},paramName);
      ConnFC.sim.spec.dimLabels = [{num2cell(1:66),num2cell(1:66),{5,11,25}},paramVar];
      ConnFC.sim.spec.dimSize = size(plv);
      ConnFC.sim.spec.metrics = {'icoh','coh','plv','cohy','lpc'};
      
      save([savepath '/ConnFC.mat'],'-v7.3','-struct','ConnFC')
      
    end
    
    function compareSimExp = compareWithEEG(this,sim)
      p = this.params.ConnectomeEnvelopeReduce;
      savepath = fullfile(this.workpath,p.outDirectory);
      if ~exist(savepath,'dir')
        mkdir(savepath);
      end
      
      %% load simulation values:
      if nargin<2
        data = load([savepath '/ConnFC.mat']);
      else
        data = sim;
      end
      paths = dataPaths( );
      
      %% load eeg values:
      if p.eegDatabase==1 || p.eegDatabase==2 || p.eegDatabase==3
        
        if p.eegDatabase==3
          dataEegTmp = load([paths.databases '/SC_Bastian/icoh_all_lcmvhilbert_20140804.mat']);
          data.eeg.spec.metrics = {'icoh','coh','plv','cohy'};
        else
          if p.eegDatabase==2
            dataEegTmp = load([paths.databases '/SC_Bastian/icoh_all_mne2_20140730.mat']);
          else
            dataEegTmp = load([paths.databases '/SC_Bastian/icoh_all_0_20140715.mat']);
          end
          data.eeg.spec.metrics = {'icoh','coh','plv'};
        end
        
        for n=1:length(data.eeg.spec.metrics)
          % replace NaNs of subject 6 on day 2 with day 1:
          dataEegTmp.([data.eeg.spec.metrics{n} '_all'])(6,2,:,:,:,:) = dataEegTmp.([data.eeg.spec.metrics{n} '_all'])(6,1,:,:,:,:);
          
           % permute to format [roi, roi, freq, more-dims...]:
          data.eeg.(data.eeg.spec.metrics{n}) = permute(dataEegTmp.([data.eeg.spec.metrics{n} '_all']),[5 6 7 1 2 3 4]);
        end
        data.eeg.spec.dimName = {'roi','roi','freq','subjEeg','day','prepost','task'};
        data.eeg.spec.dimLabels = {num2cell(1:66),num2cell(1:66),{5,11,25},num2cell(1:10),{1,2},{'pre','post'},{'co','ce'}};
        
      elseif p.eegDatabase==4 || ...
          p.eegDatabase==5 || ...
          p.eegDatabase==6 || ...
          p.eegDatabase==7 || ...
          p.eegDatabase==10 || ...
          p.eegDatabase==11 || ...
          p.eegDatabase==12 || ...
          p.eegDatabase==13 || ...
          p.eegDatabase==14 || ...
          p.eegDatabase==15 || ...
          p.eegDatabase==16 || ...
          p.eegDatabase==17 || ...
          p.eegDatabase==18 || ...
          p.eegDatabase==19 || ...
          p.eegDatabase==20 || ...
          p.eegDatabase==21 || ...
          p.eegDatabase==22 || ...
          p.eegDatabase==23 || ...
          p.eegDatabase==24 || ...
          p.eegDatabase==25 || ...
          p.eegDatabase==26 || ...
          p.eegDatabase==27 || ...
          p.eegDatabase==28 || ...
          p.eegDatabase==29
        
        if p.eegDatabase==4
          dataEegTmp = load([paths.databases '/SC_Bastian/icoh_all_lcmvhilbertrest_20140807.mat']);
        elseif p.eegDatabase==5
          dataEegTmp = load([paths.databases '/SC_Bastian/icoh_all_lcmv0_001hilbert_rest_bp_20140902.mat']);
        elseif p.eegDatabase==6
          dataEegTmp = load([paths.databases '/SC_Bastian/icoh_all_lcmv0_001_bp_hilbert_20140907.mat']);
        elseif p.eegDatabase==7
          dataEegTmp = load([paths.databases '/SC_Bastian/icoh_eeg_funcconn_bp_lcmv0_001_hilbert_20140912.mat']);
        elseif p.eegDatabase==10
          dataEegTmp = load([paths.databases '/SC_Bastian/icoh_eeg_funcconn_bp_lcmv0_001_hilbert_20141002.mat']);
        elseif p.eegDatabase==11
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141002_funcconn_bp_elor0_001_hilbert.mat']);
        elseif p.eegDatabase==12
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141003_funcconn_bp_mne_prewhit_0_001_hilbert.mat']);
        elseif p.eegDatabase==13
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141003_funcconn_bp_mne_noprewhit_0_001_hilbert.mat']);
        elseif p.eegDatabase==14
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141008_funcconn_bp_lcmv0_001_hilbert_LPC.mat']);
        elseif p.eegDatabase==15
          dataEegTmp = load([paths.databases '/SC_Bastian/20141006_conn_bp.mat']);
        elseif p.eegDatabase==16
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141028_funconn_bp11_lcmv_hilbert.mat']);
        elseif p.eegDatabase==17
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141103_controls_fs_funconn_bp_lcmv_hilbert_3_30.mat']);
        elseif p.eegDatabase==18
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141109_controls_fs_funconn_bp_pinvonall_hilbert_3_30.mat']);
        elseif p.eegDatabase==19
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141109_controls_fs_funconn_bp_pinvperroi_hilbert_3_30.mat']);
        elseif p.eegDatabase==20
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141212_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30.mat']);
        elseif p.eegDatabase==21
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150102_controls_fs_funconn_bponetrial_elor_hilbert_51122_020115.mat']);
        elseif p.eegDatabase==22
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20141124_controls_bponetrial_elor_hilbert_51122.mat']);
        elseif p.eegDatabase==23
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150113_controls_fs_funconn_lcmv_bponetrial_hilbert_3_30.mat']);
        elseif p.eegDatabase==24
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150114_controls_fs_funconn_elor_bponetrial_hilbert_51122_120115.mat']);
        elseif p.eegDatabase==25
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150114_controls_fs_funconn_lcmv_bponetrial_hilbert_3_30.mat']);
        elseif p.eegDatabase==26
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150114_conn_bponetrial_3_30.mat']);
        elseif p.eegDatabase==27
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150125_controls_fs_funconn_bponetrial_lcmv_hilbert_3_30_entmirrored.mat']);
        elseif p.eegDatabase==28
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150206_fc_patients_mne_bp_hilbert.mat']);
        elseif p.eegDatabase==29
          dataEegTmp = load([paths.databases '/SC_Bastian/eeg_20150208_fc_patients_mne0_001_bp_hilbert.mat']);
        end
        
        if p.eegDatabase==14
          data.eeg.spec.metrics = {'icoh','coh','plv','cohy','lpc'};
        elseif p.eegDatabase==16 || p.eegDatabase==17 || p.eegDatabase==20 || p.eegDatabase==21 || p.eegDatabase==22 || p.eegDatabase==23 || p.eegDatabase==24 || p.eegDatabase==25 || p.eegDatabase==27 || p.eegDatabase==28 || p.eegDatabase==29
          data.eeg.spec.metrics = {'icoh','coh','plv','cohy','lpc','pli','wpli'};
        else
          data.eeg.spec.metrics = {'icoh','coh','plv','cohy'};
        end
        
        if p.eegDatabase==15 || p.eegDatabase==26
          numRois=63;
        else
          numRois=66;
        end
        
        for n=1:length(data.eeg.spec.metrics)
          metName = data.eeg.spec.metrics{n};
          if strcmp(metName,'lpc')
            metName = 'LPC';
          end          
          if strcmp(metName,'pli')
            metName = 'PLI';
          end          
          if strcmp(metName,'wpli')
            metName = 'WPLI';
          end          
          % replace NaNs of subject 6 on day 2 with day 1:
          dataEegTmp.([metName '_all'])(6,2,:,:,:,:) = dataEegTmp.([metName '_all'])(6,1,:,:,:,:);
          % replace NaNs of subject 7 on day 2 post_rs with day 1 post_rs:
          dataEegTmp.([metName '_all'])(7,2,6,:,:,:) = dataEegTmp.([metName '_all'])(7,1,6,:,:,:);
          
          % permute to format [roi, roi, freq, more-dims...]:
          data.eeg.(data.eeg.spec.metrics{n}) = permute(dataEegTmp.([metName '_all']),[4 5 6 1 2 3]);
        end
        data.eeg.spec.dimName = {'roi','roi','freq','subjEeg','day','cond'};
        if p.eegDatabase==17 || p.eegDatabase==18 || p.eegDatabase==19 || p.eegDatabase==20 || p.eegDatabase==21 || p.eegDatabase==23
          freqz = num2cell(1:30);
          dayz = {1};
        elseif p.eegDatabase==22
          freqz = num2cell(1:22);
          dayz = {1};
        elseif p.eegDatabase==24
          freqz = num2cell(1:22);
          dayz = {1,2};
        elseif p.eegDatabase==25 || p.eegDatabase==26 || p.eegDatabase==27 || p.eegDatabase==28 || p.eegDatabase==29
          freqz = num2cell(1:30);
          dayz = {1,2};
        else
          freqz = {5,11,25};
          dayz = {1,2};
        end

        if p.eegDatabase==20 || p.eegDatabase==21 || p.eegDatabase==23 || p.eegDatabase==24 || p.eegDatabase==25 || p.eegDatabase==26 || p.eegDatabase==27 || p.eegDatabase==28 || p.eegDatabase==29
          subj=num2cell(1:20);
        else
          subj=num2cell(1:10);
        end
        data.eeg.spec.dimLabels = {num2cell(1:numRois),num2cell(1:numRois),freqz,subj,dayz,{'gripper_ce','gripper_co','glove_ce','glove_co','pre_rs','post_rs'}};
      
      elseif p.eegDatabase==8
        
        dataEegTmp = load([paths.databases '/SC_Bastian/icoh_patients_eeg_funcconn_bp_lcmv0001_hilbert_20140912.mat']);
        data.eeg.spec.metrics = {'coh','plv'};
        for n=1:length(data.eeg.spec.metrics)
          % permute to format [roi, roi, freq, more-dims...]:
          data.eeg.(data.eeg.spec.metrics{n}) = permute(dataEegTmp.([data.eeg.spec.metrics{n} '_all']),[4 5 6 1 2 3]);
        end
        data.eeg.spec.dimName = {'roi','roi','freq','subjEeg','day','cond'};
        data.eeg.spec.dimLabels = {num2cell(1:66),num2cell(1:66),{5,11,25},num2cell(1:10),{1,2,3},{'gripper_ce','gripper_co','glove_ce','glove_co','pre_rs','post_rs'}};
        
      elseif p.eegDatabase==9
        
        dataEegTmp = load([paths.databases '/SC_Bastian/patients_eeg_funcconn_bp_lcmv0001_hilbert_20140930.mat']);
        data.eeg.spec.metrics = {'icoh','coh','plv','cohy'};
        for n=1:length(data.eeg.spec.metrics)
          % permute to format [roi, roi, freq, more-dims...]:
          data.eeg.(data.eeg.spec.metrics{n}) = permute(dataEegTmp.([data.eeg.spec.metrics{n} '_all']),[4 5 6 1 2 3]);
        end
        data.eeg.spec.dimName = {'roi','roi','freq','subjEeg','day','cond'};
        data.eeg.spec.dimLabels = {num2cell(1:66),num2cell(1:66),{5,11,25},num2cell(1:10),{1,2,3},{'gripper_ce','gripper_co','glove_ce','glove_co','pre_rs','post_rs'}};

      end
      data.eeg.spec.dimSize = size(data.eeg.plv);
      
      %% pre filter and average over the data:
      data.eeg = ConnectomeEnvelopeReduce.filterTensor(data.eeg, p.eeg);
      data.sim = ConnectomeEnvelopeReduce.filterTensor(data.sim, p.sim);
      
      numRois = size(data.eeg.plv,1);
      numFreq = size(data.eeg.plv,3);
      trigIds=find(triu(ones(numRois,numRois),1));
      
      %% convert FC-Matrix to Vector:
      data.eeg = ConnectomeEnvelopeReduce.tensorCombineDims(data.eeg, 'roi', 'trig', 'roiPair');
      data.sim = ConnectomeEnvelopeReduce.tensorCombineDims(data.sim, 'roi', 'trig', 'roiPair');
      
      
      %% Add Absolute Imaginary Coherence:
      if isfield(data.eeg,'icoh')
        data.eeg.aicoh = abs(data.eeg.icoh);
        data.sim.aicoh = abs(data.sim.icoh);
      end
      
      if isfield(data.eeg,'lpc') && isfield(data.sim,'lpc')
        data.eeg.alpc = abs(data.eeg.lpc);
        data.sim.alpc = abs(data.sim.lpc);
      end
      if isfield(data.eeg,'pli') && isfield(data.sim,'pli')
        data.eeg.apli = abs(data.eeg.pli);
        data.sim.apli = abs(data.sim.pli);
      end
      if isfield(data.eeg,'wpli') && isfield(data.sim,'wpli')
        data.eeg.awpli = abs(data.eeg.wpli);
        data.sim.awpli = abs(data.sim.wpli);
      end
      
      %% save temp data
      if p.savePreEvaluationData
        save([savepath '/data.mat'],'-v7.3','-struct','data')
      end
      
      
      %% Start Evaluation:
      simMetrics = setdiff(fieldnames(data.sim),{'spec'});
      eegMetrics = setdiff(fieldnames(data.eeg),{'spec'});
      metrics = intersect(eegMetrics, simMetrics);
      compareSimExp = struct();
      
      %% combine dimNames and Labels and Sizes:
      if sum(strcmp(data.sim.spec.dimName,'freq'))
        simDimSize = data.sim.spec.dimSize(3:end);
        simDimLabels = data.sim.spec.dimLabels(3:end);
        simDimName = data.sim.spec.dimName(3:end);
      else
        simDimSize = data.sim.spec.dimSize(2:end);
        simDimLabels = data.sim.spec.dimLabels(2:end);
        simDimName = data.sim.spec.dimName(2:end);
      end
      eegDimSize = data.eeg.spec.dimSize(3:end);
      eegDimLabels = data.eeg.spec.dimLabels(3:end);
      eegDimName = data.eeg.spec.dimName(3:end);
      
      specTotal.dimSize = [numRois numRois numFreq eegDimSize simDimSize];
      specTotal.dimLabels = [{num2cell(1:numRois),num2cell(1:numRois)} data.eeg.spec.dimLabels(2) eegDimLabels simDimLabels];
      specTotal.dimName = [{'roi', 'roi', 'freq'} eegDimName simDimName];
      
      specWithoutFreq = specTotal;
      specWithoutFreq.dimSize(3) = [];
      specWithoutFreq.dimLabels(3) = [];
      specWithoutFreq.dimName(3) = [];
      
      specWithoutRois = specTotal;
      specWithoutRois.dimSize(1:2) = [];
      specWithoutRois.dimLabels(1:2) = [];
      specWithoutRois.dimName(1:2) = [];
      
      specWithoutRoisFreq = specTotal;
      specWithoutRoisFreq.dimSize(1:3) = [];
      specWithoutRoisFreq.dimLabels(1:3) = [];
      specWithoutRoisFreq.dimName(1:3) = [];
      
      specOneRoisFreq = specTotal;
      specOneRoisFreq.dimSize(1) = [];
      specOneRoisFreq.dimLabels(1) = [];
      specOneRoisFreq.dimName(1) = [];
      
      specOneRois = specWithoutFreq;
      specOneRois.dimSize(1) = [];
      specOneRois.dimLabels(1) = [];
      specOneRois.dimName(1) = [];
      
      if p.calcFisher && p.calcFisher<3
        for m=1:length(metrics)
          if p.calcFisher==2
            if isreal(data.eeg.(metrics{m}))
              data.eeg.(metrics{m}) = atanh(data.eeg.(metrics{m}));
            end
            if isreal(data.sim.(metrics{m}))
              data.sim.(metrics{m}) = atanh(data.sim.(metrics{m}));
            end
          else
            if isreal(data.eeg.(metrics{m}))
              if sum(data.eeg.(metrics{m})(:)<=0)
                data.eeg.(metrics{m}) = atanh(data.eeg.(metrics{m}));
              else
                data.eeg.(metrics{m}) = log(data.eeg.(metrics{m}));
              end
            end
            if isreal(data.sim.(metrics{m}))
              if sum(data.sim.(metrics{m})(:)<=0)
                data.sim.(metrics{m}) = atanh(data.sim.(metrics{m}));
              else
                data.sim.(metrics{m}) = log(data.sim.(metrics{m}));
              end
            end
          end
        end
      end
      
      
      if p.calcCompareAicohWithCoh
        metrics{end+1} = 'aicohcoh';
      end
      
      numTrig = numRois*(numRois-1)/2;
      
      %% per freq:
      for m=1:length(metrics)
        
        disp(['metric : ' num2str(m) '/' num2str(length(metrics))])
        if ~strcmp(metrics{m},'cohy') && ~strcmp(metrics{m},'aicohcoh')
          if p.calcSquaredDist
            compareSimExp.distRoiPair.sqrDist.perFreq.(metrics{m}) = zeros(numRois,numRois,numFreq,prod(eegDimSize),prod(simDimSize));
            compareSimExp.distRoiPair.absDist.perFreq.(metrics{m}) = zeros(numRois,numRois,numFreq,prod(eegDimSize),prod(simDimSize));
          end
          if p.calcSquaredDistAvg
            compareSimExp.distAvgPerRoi.sqrDist.perFreq.(metrics{m}) = zeros(numRois,numFreq,prod(eegDimSize),prod(simDimSize));
            compareSimExp.distAvgPerRoi.absDist.perFreq.(metrics{m}) = zeros(numRois,numFreq,prod(eegDimSize),prod(simDimSize));
          end
        end
        
        for freq=1:numFreq
          
          if strcmp(metrics{m},'aicohcoh')
            eegMetric = 'aicoh';
            simMetric = 'coh';
          else
            eegMetric = metrics{m};
            simMetric = metrics{m};
          end
          
          
          
          
          
          tmpEEG = reshape(data.eeg.(eegMetric)(:,freq,:,:,:,:,:),[numTrig prod(eegDimSize)]);
          if sum(strcmp(data.sim.spec.dimName,'freq'))
            tmpSIM = reshape(data.sim.(simMetric)(:,freq,:,:,:,:,:),[numTrig prod(simDimSize)]);
          else
            tmpSIM = reshape(data.sim.(simMetric)(:,:,:,:,:,:),[numTrig prod(simDimSize)]);
          end
          
          
          
          if p.calcCrossFreq
            for freqSim=1:numFreq
              tmpSIM_otherFreq = reshape(data.sim.(simMetric)(:,freqSim,:,:,:,:,:),[numTrig prod(simDimSize)]);
              if strcmp(metrics{m},'cohy')
                compareSimExp.crossFreq.(metrics{m}).coh(freq,freqSim,:,:) = ConnectomeEnvelopeReduce.calcCoherence(tmpEEG',tmpSIM_otherFreq');
              else
                %% Corr:
                if p.calcFisher==3
                  [compareSimExp.crossFreq.(metrics{m}).rho(freq,freqSim,:,:), compareSimExp.crossFreq.(metrics{m}).pval(freq,freqSim,:,:)] = corr(tmpEEG,tmpSIM_otherFreq,'type','Spearman');
                else
                  [compareSimExp.crossFreq.(metrics{m}).rho(freq,freqSim,:,:), compareSimExp.crossFreq.(metrics{m}).pval(freq,freqSim,:,:)] = corr(tmpEEG,tmpSIM_otherFreq);
                end
              end
            end
          end
          
          if strcmp(metrics{m},'cohy')
            compareSimExp.perFreq.(metrics{m}).coh(freq,:,:) = ConnectomeEnvelopeReduce.calcCoherence(tmpEEG',tmpSIM');
          else
            %% Corr:
            if p.calcFisher==3
              [compareSimExp.perFreq.(metrics{m}).rho(freq,:,:), compareSimExp.perFreq.(metrics{m}).pval(freq,:,:)] = corr(tmpEEG,tmpSIM,'type','Spearman');
            else
              [compareSimExp.perFreq.(metrics{m}).rho(freq,:,:), compareSimExp.perFreq.(metrics{m}).pval(freq,:,:)] = corr(tmpEEG,tmpSIM);
            end
            
            %% Partial Corr:
            if p.calcPartialCorrEuclDist
              paths = dataPaths();
              dist = load(fullfile(paths.databases,'SC_Bastian','icoh_all_lcmvhilbertrest_20140807.mat'),'pos');
              dist = squeeze(sqrt(sum(bsxfun(@minus,dist.pos,permute(dist.pos,[3 2 1])).^2,2)));
              dist = dist(trigIds);
              [compareSimExp.perFreq.(metrics{m}).partialEuclDistRho(freq,:,:), compareSimExp.perFreq.(metrics{m}).partialEuclDistPval(freq,:,:)] = partialcorr(tmpEEG,tmpSIM,dist);
            end
            if p.calcPartialCorrFiberDist
              paths = dataPaths();
              dist = load(fullfile(paths.databases,'SC_Bastian','dist_and_CI_controls_preprocessed.mat'),'avg_dist');
              dist = dist.avg_dist(trigIds);
              [compareSimExp.perFreq.(metrics{m}).partialFiberDistRho(freq,:,:), compareSimExp.perFreq.(metrics{m}).partialFiberDistPval(freq,:,:)] = partialcorr(tmpEEG,tmpSIM,dist);
            end
            
            %% Jaccard similarity coefficient (only for positive real valued metrics):
            if ~strcmp(metrics{m},'icoh')
              tmpMin = bsxfun(@min,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
              tmpMax = bsxfun(@max,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
              compareSimExp.perFreq.(metrics{m}).jac(freq,:,:) = sum(tmpMin,1) ./ sum(tmpMax,1);
            end
            
            %% Linear Regression using Total Least-Squared:
            if p.calcSquaredDist || p.calcSquaredDistAvg
              for eegId=1:size(tmpEEG,2)
                for simId=1:size(tmpSIM,2)
                  [~, ~, d] = fit_2D_data(tmpSIM(:,simId), tmpEEG(:,eegId), 'no');
                  squaredDist = zeros(numRois,numRois);
                  squaredDist(trigIds) = d;
                  
                  absDist = zeros(numRois,numRois);
                  absDist(trigIds) = abs(tmpSIM(:,simId) - tmpEEG(:,eegId));
                  
                  if p.calcSquaredDist
                    compareSimExp.distRoiPair.sqrDist.perFreq.(metrics{m})(:,:,freq,eegId,simId) = squaredDist;
                    compareSimExp.distRoiPair.absDist.perFreq.(metrics{m})(:,:,freq,eegId,simId) = absDist;
                  end
                  if p.calcSquaredDistAvg
                    compareSimExp.distAvgPerRoi.sqrDist.perFreq.(metrics{m})(:,freq,eegId,simId) = mean(squaredDist(:,2:end) + squaredDist(1:end-1,:)',2);
                    compareSimExp.distAvgPerRoi.absDist.perFreq.(metrics{m})(:,freq,eegId,simId) = mean(absDist(:,2:end) + absDist(1:end-1,:)',2);
                  end
                  
                  
                  
                end
              end
              
            end
            
            
            
          end
        end
        
        fnames = fieldnames(compareSimExp.perFreq.(metrics{m}));
        for f=1:length(fnames)
          compareSimExp.perFreq.(metrics{m}).(fnames{f}) = reshape(compareSimExp.perFreq.(metrics{m}).(fnames{f}), [numFreq eegDimSize simDimSize]);
          compareSimExp.perFreqAvg.(metrics{m}).(fnames{f}) = permute(mean(compareSimExp.perFreq.(metrics{m}).(fnames{f}),1),[2:length(size(compareSimExp.perFreq.(metrics{m}).(fnames{f}))) 1]);
        end
        
        
        if p.calcCrossFreq
          fnames = fieldnames(compareSimExp.crossFreq.(metrics{m}));
          for f=1:length(fnames)
            compareSimExp.crossFreq.(metrics{m}).(fnames{f}) = reshape(compareSimExp.crossFreq.(metrics{m}).(fnames{f}), [numFreq numFreq eegDimSize simDimSize]);
          end
        end
        
        if ~strcmp(metrics{m},'cohy')
          if p.calcSquaredDist
            compareSimExp.distRoiPair.sqrDist.perFreq.(metrics{m}) = reshape(compareSimExp.distRoiPair.sqrDist.perFreq.(metrics{m}),[numRois numRois numFreq eegDimSize simDimSize]);
            compareSimExp.distRoiPair.sqrDist.perFreqAvg.(metrics{m}) = permute(mean(compareSimExp.distRoiPair.sqrDist.perFreq.(metrics{m}),3),[1 2 4:length(size(compareSimExp.distRoiPair.sqrDist.perFreq.(metrics{m}))) 3]);
            compareSimExp.distRoiPair.absDist.perFreq.(metrics{m}) = reshape(compareSimExp.distRoiPair.absDist.perFreq.(metrics{m}),[numRois numRois numFreq eegDimSize simDimSize]);
            compareSimExp.distRoiPair.absDist.perFreqAvg.(metrics{m}) = permute(mean(compareSimExp.distRoiPair.absDist.perFreq.(metrics{m}),3),[1 2 4:length(size(compareSimExp.distRoiPair.absDist.perFreq.(metrics{m}))) 3]);
          end
          if p.calcSquaredDistAvg
            compareSimExp.distAvgPerRoi.sqrDist.perFreq.(metrics{m}) = reshape(compareSimExp.distAvgPerRoi.sqrDist.perFreq.(metrics{m}),[numRois numFreq eegDimSize simDimSize]);
            compareSimExp.distAvgPerRoi.sqrDist.perFreqAvg.(metrics{m}) = permute(mean(compareSimExp.distAvgPerRoi.sqrDist.perFreq.(metrics{m}),2),[1 3:length(size(compareSimExp.distAvgPerRoi.sqrDist.perFreq.(metrics{m}))) 2]);
            compareSimExp.distAvgPerRoi.absDist.perFreq.(metrics{m}) = reshape(compareSimExp.distAvgPerRoi.absDist.perFreq.(metrics{m}),[numRois numFreq eegDimSize simDimSize]);
            compareSimExp.distAvgPerRoi.absDist.perFreqAvg.(metrics{m}) = permute(mean(compareSimExp.distAvgPerRoi.absDist.perFreq.(metrics{m}),2),[1 3:length(size(compareSimExp.distAvgPerRoi.absDist.perFreq.(metrics{m}))) 2]);
          end
        end
      end
      
      %% add specs:
      compareSimExp.perFreq.spec = specWithoutRois;
      compareSimExp.perFreqAvg.spec = specWithoutRoisFreq;
      if p.calcSquaredDist
        compareSimExp.distRoiPair.sqrDist.perFreq.spec = specTotal;
        compareSimExp.distRoiPair.sqrDist.perFreqAvg.spec = specWithoutFreq;
        compareSimExp.distRoiPair.absDist.perFreq.spec = specTotal;
        compareSimExp.distRoiPair.absDist.perFreqAvg.spec = specWithoutFreq;
      end
      if p.calcSquaredDistAvg
        compareSimExp.distAvgPerRoi.sqrDist.perFreq.spec = specOneRoisFreq;
        compareSimExp.distAvgPerRoi.sqrDist.perFreqAvg.spec = specOneRois;
        compareSimExp.distAvgPerRoi.absDist.perFreq.spec = specOneRoisFreq;
        compareSimExp.distAvgPerRoi.absDist.perFreqAvg.spec = specOneRois;
      end
      if p.calcCrossFreq
        compareSimExp.crossFreq.spec = specWithoutRois;
        compareSimExp.crossFreq.spec.dimSize = [3 compareSimExp.crossFreq.spec.dimSize];
        compareSimExp.crossFreq.spec.dimLabels = [compareSimExp.crossFreq.spec.dimLabels(1) compareSimExp.crossFreq.spec.dimLabels];
        compareSimExp.crossFreq.spec.dimName = [{'freq EEG','freq SIM'} compareSimExp.crossFreq.spec.dimName(2:end)];
      end
      
      eegSimDimSize = [eegDimSize simDimSize];
      if length(eegSimDimSize)<2
        eegSimDimSize(2) = 1;
      end
      
      %% combine all freq:
      for m=1:length(metrics)
        
        if strcmp(metrics{m},'aicohcoh')
          eegMetric = 'aicoh';
          simMetric = 'coh';
        else
          eegMetric = metrics{m};
          simMetric = metrics{m};
        end
        
        if sum(strcmp(data.sim.spec.dimName,'freq'))
          tmpEEG = reshape(data.eeg.(eegMetric),[numTrig*numFreq prod(eegDimSize)]);
          tmpSIM = reshape(data.sim.(simMetric),[numTrig*numFreq prod(simDimSize)]);
        else
          tmpEEG = reshape(mean(data.eeg.(eegMetric),2),[numTrig prod(eegDimSize)]);
          tmpSIM = reshape(data.sim.(simMetric),[numTrig prod(simDimSize)]);
        end
        
        if strcmp(metrics{m},'cohy')
          compareSimExp.overFreq.(metrics{m}).coh = ConnectomeEnvelopeReduce.calcCoherence(tmpEEG',tmpSIM');
          compareSimExp.overFreq.(metrics{m}).coh = reshape(compareSimExp.overFreq.(metrics{m}).coh,eegSimDimSize);
        else
          %% Corr:
          [compareSimExp.overFreq.(metrics{m}).rho, compareSimExp.overFreq.(metrics{m}).pval] = corr(tmpEEG,tmpSIM);
          compareSimExp.overFreq.(metrics{m}).rho = reshape(compareSimExp.overFreq.(metrics{m}).rho,eegSimDimSize);
          compareSimExp.overFreq.(metrics{m}).pval = reshape(compareSimExp.overFreq.(metrics{m}).pval,eegSimDimSize);
          
          %% Jaccard similarity coefficient (only for positive real valued metrics):
          if ~strcmp(metrics{m},'icoh')
            tmpMin = bsxfun(@min,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
            tmpMax = bsxfun(@max,permute(tmpEEG,[1 2]),permute(tmpSIM,[1 3 2]));
            compareSimExp.overFreq.(metrics{m}).jac = sum(tmpMin,1) ./ sum(tmpMax,1);
            compareSimExp.overFreq.(metrics{m}).jac = reshape(compareSimExp.overFreq.(metrics{m}).jac, eegSimDimSize);
          end
        end
      end
      compareSimExp.overFreq.spec = specWithoutRoisFreq;

      
      %% maybe reduce some dimensions in the results:
      if isfield(p,'results')
        compareSimExp.overFreq = ConnectomeEnvelopeReduce.filterTensor(compareSimExp.overFreq, p.results);
        compareSimExp.perFreqAvg = ConnectomeEnvelopeReduce.filterTensor(compareSimExp.perFreqAvg, p.results);
        if p.calcSquaredDist
          compareSimExp.distRoiPair = ConnectomeEnvelopeReduce.filterTensor(compareSimExp.distRoiPair, p.results);
        end
        if p.calcSquaredDistAvg
          compareSimExp.distAvgPerRoi = ConnectomeEnvelopeReduce.filterTensor(compareSimExp.distAvgPerRoi, p.results);
        end
      end
      
      
      %% maybe normalize the paired subject matrizes:
      if p.doNormalizeSubjectMat
        compareSimExp.overFreq = ConnectomeEnvelopeReduce.normalizeMat(compareSimExp.overFreq, [], '' );
        compareSimExp.perFreqAvg = ConnectomeEnvelopeReduce.normalizeMat(compareSimExp.perFreqAvg, [], '' );        
      end
      
      
      %% evaluate subject specific permutation tests:
      if p.calcSubjectSpecificTests && sum(strcmp(compareSimExp.overFreq.spec.dimName,'subjEeg')) && sum(strcmp(compareSimExp.overFreq.spec.dimName,'subjId'))
        compareSimExp.permTest.overFreq = ConnectomeEnvelopeReduce.evalPermTest(compareSimExp.overFreq, [], '^pval' );
        compareSimExp.permTest.perFreqAvg = ConnectomeEnvelopeReduce.evalPermTest(compareSimExp.perFreqAvg, [], '^pval' );
      end
      
      
      
      
      %% evaluate phase shifts in sim and eeg:
      if p.calcPhaseLags && sum(strcmp(metrics,'cohy')) && sum(strcmp(data.sim.spec.dimName,'freq'))
        
        tmpEEG = reshape(data.eeg.cohy,[numTrig numFreq prod(eegDimSize)]);
        tmpSIM = reshape(data.sim.cohy,[numTrig numFreq prod(simDimSize)]);
        
        phasediffEEG = angle(tmpEEG);
        phasediffSIM = angle(tmpSIM);
        
        numBins = 15;
        
        phaseBinEdges = linspace(-pi,pi,numBins+1);
        phaseBins = phaseBinEdges(1:end-1)+diff(phaseBinEdges);
        [nEEG,binEEG] = histc(phasediffEEG,phaseBinEdges,1);
        [nSIM,binSim] = histc(phasediffSIM,phaseBinEdges,1);
        
        nEEG(end-1,:,:) = nEEG(end-1,:,:) + nEEG(end,:,:);
        nSIM(end-1,:,:) = nSIM(end-1,:,:) + nSIM(end,:,:);
        
        nEEG(end,:,:) = [];
        nSIM(end,:,:) = [];
        
        compareSimExp.phaseBins.eeg.nEEG = reshape(nEEG,[size(nEEG,1) size(nEEG,2) eegDimSize]);
        compareSimExp.phaseBins.eeg.spec.dimSize = size(compareSimExp.phaseBins.eeg.nEEG);
        compareSimExp.phaseBins.eeg.spec.dimLabels = [{num2cell(phaseBins)} data.eeg.spec.dimLabels(2) eegDimLabels];
        compareSimExp.phaseBins.eeg.spec.dimName = [{'phase lags', 'freq'} eegDimName];
        
        compareSimExp.phaseBins.sim.nSIM = reshape(nSIM,[size(nEEG,1) size(nEEG,2) simDimSize]);
        compareSimExp.phaseBins.sim.spec.dimSize = size(compareSimExp.phaseBins.sim.nSIM);
        compareSimExp.phaseBins.sim.spec.dimLabels = [{num2cell(phaseBins)} data.sim.spec.dimLabels(2) simDimLabels];
        compareSimExp.phaseBins.sim.spec.dimName = [{'phase lags', 'freq'} simDimName];        
        
        histSp = bsxfun(@times,permute(nEEG,[1 4 2 3 5]),permute(nSIM,[4 1 2 5 3]));
        hist2D = zeros(numBins,numBins,3,prod(eegDimSize),prod(simDimSize));
        for eegId=1:size(tmpEEG,3)
          for simId=1:size(tmpSIM,3)
            for freq=1:3
              hist2D(:,:,freq,eegId,simId) = accumarray([binEEG(:,freq,eegId) binSim(:,freq,simId)],ones(size(binSim,1),1));
            end
          end
        end
        
        compareSimExp.phaseLags2DHist.histSp = bsxfun(@rdivide,histSp,sum(sum(histSp,1),2));
        compareSimExp.phaseLags2DHist.hist2D = bsxfun(@rdivide,hist2D,sum(sum(hist2D,1),2));
        
        compareSimExp.phaseLags2DHist.spec = specTotal;
        compareSimExp.phaseLags2DHist.spec.dimSize(1:2) = [numBins numBins];
        compareSimExp.phaseLags2DHist.spec.dimLabels(1:2) = {num2cell(phaseBins), num2cell(phaseBins)};
        compareSimExp.phaseLags2DHist.spec.dimName(1:2) = {'eegPhaseLag', 'simPhaseLag'};
        
        tmp = compareSimExp.phaseLags2DHist.hist2D .* log2( compareSimExp.phaseLags2DHist.hist2D ./ compareSimExp.phaseLags2DHist.histSp );
        tmp(isnan(tmp)) = 0;
        compareSimExp.phaseLagsEval.KL = sum(sum( tmp ,1),2);        
        compareSimExp.phaseLagsEval.KL = reshape(compareSimExp.phaseLagsEval.KL,[3 eegDimSize simDimSize]);

        compareSimExp.phaseLagsEval.spec = specWithoutRois;
        
        %% calc the following only for the simulation with the best coh value:
        [~,Ibest] = max(compareSimExp.perFreqAvg.coh.rho(:));
        if isempty(eegDimSize)
          compareSimExp.samples.bestEegId = 1;
          compareSimExp.samples.bestSimId = Ibest;
        elseif isempty(simDimSize)
          compareSimExp.samples.bestEegId = Ibest;
          compareSimExp.samples.bestSimId = 1;
        else
          [compareSimExp.samples.bestEegId,compareSimExp.samples.bestSimId] = ind2sub(eegSimDimSize,Ibest);
        end
        [I(1),I(2),I(3),I(4),I(5)] = ind2sub(compareSimExp.perFreqAvg.spec.dimSize,Ibest);
        compareSimExp.samples.I = I;
        
        compareSimExp.samples.Ibest = Ibest;
        compareSimExp.samples.cohy_eegBest = tmpEEG(:,:,compareSimExp.samples.bestEegId);
        compareSimExp.samples.cohy_simBest = tmpSIM(:,:,compareSimExp.samples.bestSimId);
        
        compareSimExp.samples.histAngles = linspace(-pi,pi,100);
        compareSimExp.samples.histAngleCohy = histc(angle(tmpEEG(:)),compareSimExp.samples.histAngles)/size(tmpEEG,1);
        
        
        %% calc desc:
        dimNames = compareSimExp.perFreqAvg.spec.dimName;
        dimLabels = compareSimExp.perFreqAvg.spec.dimLabels;
        compareSimExp.samples.descParam = '';
        for ii=1:min(length(dimNames),length(I))
          compareSimExp.samples.descParam = [compareSimExp.samples.descParam dimNames{ii} '=' num2str(dimLabels{ii}{I(ii)}) ' '];
        end
        
      end
      
      
      %% maybe combine subject dimensions from eeg and sim into one dimension:
      if strcmp(p.resultsCombineSubjDims,'match') || strcmp(p.resultsCombineSubjDims,'nonmatch') || strcmp(p.resultsCombineSubjDims,'avgMatch')
        compareSimExp.overFreq = ConnectomeEnvelopeReduce.tensorCombineDims(compareSimExp.overFreq, {'subjEeg','subjId'}, p.resultsCombineSubjDims, 'subject');
        compareSimExp.perFreqAvg = ConnectomeEnvelopeReduce.tensorCombineDims(compareSimExp.perFreqAvg, {'subjEeg','subjId'}, p.resultsCombineSubjDims, 'subject');
        if p.calcSquaredDist
          compareSimExp.distRoiPair = ConnectomeEnvelopeReduce.tensorCombineDims(compareSimExp.distRoiPair, {'subjEeg','subjId'}, p.resultsCombineSubjDims, 'subject');
        end
        if p.calcSquaredDistAvg
          compareSimExp.distAvgPerRoi = ConnectomeEnvelopeReduce.tensorCombineDims(compareSimExp.distAvgPerRoi, {'subjEeg','subjId'}, p.resultsCombineSubjDims, 'subject');
        end
        if p.calcPhaseLags
          compareSimExp.phaseLags2DHist = ConnectomeEnvelopeReduce.tensorCombineDims(compareSimExp.phaseLags2DHist, {'subjEeg','subjId'}, p.resultsCombineSubjDims, 'subject');
          compareSimExp.phaseLagsEval = ConnectomeEnvelopeReduce.tensorCombineDims(compareSimExp.phaseLagsEval, {'subjEeg','subjId'}, p.resultsCombineSubjDims, 'subject');
        end
      end
      
      
      
      %% save everything:
      compareSimExp.sim.spec.dimName = data.sim.spec.dimName;
      compareSimExp.sim.spec.dimLabels = data.sim.spec.dimLabels;
      if isfield(data.sim,'meanOrderParam')
        compareSimExp.sim.spec.metrics.orderParamMean = data.sim.meanOrderParam;
        compareSimExp.sim.spec.metrics.orderParamStd = data.sim.stdOrderParam;
      end
      save(fullfile(savepath,['compareSimExp' num2str(this.currJobid)]),'-v7.3','-struct','compareSimExp')
      
      
      
    end
    
    
    
    function plotResults(this,compareSimExp)
      
      p = this.params.ConnectomeEnvelopeReduce;
      savepath = fullfile(this.workpath,p.outDirectory);
      if isempty(this.params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder)
        plotDir = fullfile(savepath,['plots' num2str(this.currJobid)]);
      else
        plotDir = fullfile(savepath,this.params.ConnectomeEnvelopeReduce.outDirectoryPlotFolder);
      end
      mkdir(plotDir)
      
      if nargin<2
        compareSimExp = load(fullfile(savepath,['compareSimExp' num2str(this.currJobid)]));
      end
      
      
      %% plot normal result metrics:
      ConnectomeEnvelopeReduce.plotStructure(compareSimExp.overFreq, compareSimExp.overFreq.spec, 'result.overFreq', 'result_overFreq', plotDir, {'pval','partialFiberDistPval','partialEuclDistPval'}, p.permutePlotDims, p.excludePlotFieldsRegexp)
      ConnectomeEnvelopeReduce.plotStructure(compareSimExp.perFreqAvg, compareSimExp.perFreqAvg.spec, 'result.perFreqAvg', 'result_perFreqAvg', plotDir, {'pval','partialFiberDistPval','partialEuclDistPval'}, p.permutePlotDims, p.excludePlotFieldsRegexp)
      
      ConnectomeEnvelopeReduce.plotStructure(compareSimExp.perFreq, compareSimExp.perFreq.spec, 'result.perFreq', 'result_perFreq', plotDir, {'pval','partialFiberDistPval','partialEuclDistPval'}, p.permutePlotDimsPerFreq, p.excludePlotFieldsRegexp)
      
      %% plot over fixed subject if two remaining dimensions correspond to eeg and mri subject
      if length(compareSimExp.overFreq.spec.dimName)==2 && sum(strcmp(compareSimExp.overFreq.spec.dimName,'subjEeg')) && sum(strcmp(compareSimExp.overFreq.spec.dimName,'subjId'))
        ConnectomeEnvelopeReduce.plotStructureSubj(compareSimExp.overFreq, compareSimExp.overFreq.spec, 'fixSubj.overFreq', 'fixSubj_overFreq', plotDir, {'pval','partialFiberDistPval','partialEuclDistPval'})
        ConnectomeEnvelopeReduce.plotStructureSubj(compareSimExp.perFreqAvg, compareSimExp.perFreqAvg.spec, 'fixSubj.perFreqAvg', 'fixSubj_perFreqAvg', plotDir, {'pval','partialFiberDistPval','partialEuclDistPval'})
        
        
        %% TODO: use dummyvar and glmfit...
                  % x = dummyvar([subjEEG subjSIM]);
                  % glmfit()
                  
      end
      
      %% plot best params if there is at least one remaining subject dimension:
      if p.plotBestParamPerSubj && sum(strcmp(compareSimExp.overFreq.spec.dimName,'subjEeg')) || sum(strcmp(compareSimExp.overFreq.spec.dimName,'subjId'))
%         ConnectomeEnvelopeReduce.plotBestParams(compareSimExp.overFreq, 'bestParams.overFreq', 'bestParams_overFreq', plotDir, {'pval'});
%         ConnectomeEnvelopeReduce.plotBestParams(compareSimExp.perFreqAvg, 'bestParams.overFreq', 'bestParams_overFreq', plotDir, {'pval'});
      end
      
      %% plot permutation eval:
      if p.plotPermTests && isfield(compareSimExp,'permTest')
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.permTest.overFreq, [], 'ttest.overFreq', 'ttest_overFreq', plotDir, {'anovan_pval_subjDTI','anovan_pval_subjEEG'}, p.permutePlotDimsPermTest, p.excludePlotFieldsRegexp);
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.permTest.perFreqAvg, [], 'ttest.perFreqAvg', 'ttest_perFreqAvg', plotDir, {'anovan_pval_subjDTI','anovan_pval_subjEEG'}, p.permutePlotDimsPermTest, p.excludePlotFieldsRegexp);
      end
      
      %% write permutation eval to file:
      if isfield(compareSimExp,'permTest') && isempty(compareSimExp.permTest.perFreqAvg.spec.dimName)
        fileID = fopen(fullfile(plotDir,'subjPermTest_pval.txt'),'w');
        ConnectomeEnvelopeReduce.dispStructure(compareSimExp.permTest.overFreq, 'permTest.overFreq', fileID, {'anovan_pval_subjDTI','anovan_pval_subjEEG'});
        ConnectomeEnvelopeReduce.dispStructure(compareSimExp.permTest.perFreqAvg, 'permTest.perFreqAvg', fileID, {'anovan_pval_subjDTI','anovan_pval_subjEEG'});
        fclose(fileID);
      end
      
      %% squared distances
      paths = dataPaths();
      dataSC = load([paths.databases '/SC_Bastian/dist_and_CI_controls_preprocessed.mat']);
      
      if p.calcSquaredDist        
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.distRoiPair, [], 'distRoiPair', 'distRoiPair', plotDir, {''}, p.permutePlotDimsSquaredDist, p.excludePlotFieldsRegexp);
        % compare with sc_avg_dist:
        ConnectomeEnvelopeReduce.plotCompareDistWithSC(compareSimExp.distRoiPair,[],dataSC,'compareSC.distRoiPair', 'compareSC_distRoiPair', plotDir, p.excludePlotFieldsRegexp)
      end
      
      if p.calcSquaredDistAvg
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.distAvgPerRoi, [], 'distAvgPerRoi', 'distAvgPerRoi', plotDir, {''}, p.permutePlotDimsSquaredDistAvg, p.excludePlotFieldsRegexp);
        % compare with sc_avg_roi_size:
        ConnectomeEnvelopeReduce.plotCompareDistWithSC(compareSimExp.distAvgPerRoi,[],dataSC,'compareSC.distAvgPerRoi', 'compareSC_distAvgPerRoi', plotDir, p.excludePlotFieldsRegexp)
      end
      
      
      %%
      if p.calcCrossFreq
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.crossFreq, compareSimExp.crossFreq.spec, 'result.crossFreq', 'result_crossFreq', plotDir, {'pval','partialFiberDistPval','partialEuclDistPval'}, p.permutePlotDimsCrossFreq, p.excludePlotFieldsRegexp);
      end
      
      %% plot phaselag hist
      if p.calcPhaseLags
        if length(compareSimExp.phaseLags2DHist.spec.dimSize)<=4
          ConnectomeEnvelopeReduce.plotStructure(compareSimExp.phaseLags2DHist, compareSimExp.phaseLags2DHist.spec, 'result.phaseLags2DHist', 'result_phaseLags2DHist', plotDir, {}, p.permutePlotDimsPhaseLags2DHist, p.excludePlotFieldsRegexp);
        else
            
          tmp = compareSimExp.phaseLags2DHist;
          tmp.histSp = reshape(tmp.histSp,[tmp.spec.dimSize(1:3) prod(tmp.spec.dimSize(4:end))]);
          tmp.histSp = tmp.histSp(:,:,:,compareSimExp.samples.Ibest);
          tmp.hist2D = reshape(tmp.hist2D,[tmp.spec.dimSize(1:3) prod(tmp.spec.dimSize(4:end))]);
          tmp.hist2D = tmp.hist2D(:,:,:,compareSimExp.samples.Ibest);
          
          tmp.spec.dimSize(4:end) = [];
          tmp.spec.dimLabels(4:end) = [];
          tmp.spec.dimName(4:end) = [];
          ConnectomeEnvelopeReduce.plotStructure(tmp, tmp.spec, 'result.phaseLags2DHist', 'result_phaseLags2DHist', plotDir, {}, p.permutePlotDimsPhaseLags2DHist, p.excludePlotFieldsRegexp);
        end
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.phaseLagsEval, compareSimExp.phaseLagsEval.spec, 'result.phaseLagsEval', 'result_phaseLagsEval', plotDir, {}, p.permutePlotDimsPhaseLagsEval, p.excludePlotFieldsRegexp);
        
        
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.phaseBins.eeg, compareSimExp.phaseBins.eeg.spec, 'result.phaseBins.eeg', 'result_phaseBins_eeg', plotDir, {}, p.permutePlotDimsPhaseBinsEeg, p.excludePlotFieldsRegexp);
        ConnectomeEnvelopeReduce.plotStructure(compareSimExp.phaseBins.sim, compareSimExp.phaseBins.sim.spec, 'result.phaseBins.sim', 'result_phaseBins_sim', plotDir, {}, p.permutePlotDimsPhaseBinsSim, p.excludePlotFieldsRegexp);
      end
      
      
      
      
      %% plot coh of cohy
      if isfield(compareSimExp,'samples') && isfield(compareSimExp.samples,'cohy_eegBest') && isfield(compareSimExp.samples,'cohy_simBest')
        cohy_eegBest = compareSimExp.samples.cohy_eegBest;
        cohy_simBest = compareSimExp.samples.cohy_simBest;
        
        descParam = compareSimExp.samples.descParam;
        
        %% plots:
        sfigure(1); clf;
        plotFilename = 'cohy_phasediff_bestSim';
        bothCoh=abs(cohy_simBest).*abs(cohy_eegBest);
        phaseDiffDirect=mod(angle(cohy_simBest)-angle(cohy_eegBest)+pi,2*pi)-pi;
        plot(bothCoh,phaseDiffDirect,'.','MarkerSize', 12)
        title(['Phasediff vs coh for ' descParam])
        xlabel('eegCoh * simCoh')
        ylabel('eeg-sim-difference in phase-lag between rois [rad]')
        legend('5 Hz','11 Hz','25 Hz')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');
        
        sfigure(1); clf;
        plotFilename = 'cohy_absPhasediff_bestSim';
        bothCoh=abs(cohy_simBest).*abs(cohy_eegBest);
        phaseDiffDirect=mod(angle(cohy_simBest)-angle(cohy_eegBest)+pi,2*pi)-pi;
        plot(bothCoh,abs(phaseDiffDirect),'.','MarkerSize', 12)
        title(['Phasediff vs coh for ' descParam])
        xlabel('eegCoh * simCoh')
        ylabel('abs eeg-sim-difference in phase-lag between rois [rad]')
        legend('5 Hz','11 Hz','25 Hz')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');
        
        %%
%         sfigure(1); clf;
%         plotFilename = 'hist_phaseLag_EEG_weighted';
%         bothCoh=abs(cohy_eegBest);
%         phaseDiffDirect=mod(angle(cohy_eegBest)+pi,2*pi)-pi;
%         edges=linspace(-pi,pi,30);
%         [trash bin]=histc(phaseDiffDirect,edges);
%         count=accumarray(bin,bothCoh(:));
%         bar(edges(1:end-1)+diff(edges)/2,count)
%         set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
%         export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');
% %%
%         sfigure(1); clf;
%         plotFilename = 'hist_phaseLag_SIM_weighted';
%         bothCoh=abs(cohy_simBest);
%         phaseDiffDirect=mod(angle(cohy_simBest)+pi,2*pi)-pi;
%         edges=linspace(-pi,pi,30);
%         [trash bin]=histc(phaseDiffDirect,edges);
%         count=accumarray(bin,bothCoh(:));
%         bar(edges(1:end-1)+diff(edges)/2,count)
%         set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
%         export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');

        %%
        sfigure(1); clf;
        plotFilename = ['hist_phaseLag_EEG_allfreq'];
        hist(angle(cohy_eegBest(:)),100)
        title(['Histogram over all ROI-pairs eeg'])
        xlabel('phase lag [rad]')
        ylabel('number of ROI-pairs')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');
        
        sfigure(1); clf;
        plotFilename = ['hist_phaseLag_bestSim_allfreq'];
        hist(angle(cohy_simBest(:)),100)
        title(['Histogram over all ROI-pairs for sim'])
        xlabel('phase lag [rad]')
        ylabel('number of ROI-pairs')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');
        
        
        %%
       
        edges = linspace(-pi, pi, 30);
        
        sfigure(1); clf;
        plotFilename = ['hist_phaseLag_EEG_freq'];
        [n,bin] = histc(angle(cohy_eegBest),edges,1);
        n(end,:)=[];
        %         bar(compareSimExp.samples.histAngles,compareSimExp.samples.histAngleCohy)
        plot(edges(1:end-1)+diff(edges),n)
        title(['Histogram over all ROI-pairs eeg'])
        xlabel('phase lag [rad]')
        ylabel('number of ROI-pairs')
        legend('5 Hz','11 Hz','25 Hz')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');
        
        sfigure(1); clf;
        plotFilename = ['hist_phaseLag_bestSim'];
        [n,bin] = histc(angle(cohy_simBest),edges,1);
        n(end,:)=[];
        plot(edges(1:end-1)+diff(edges),n)
        title(['Histogram over all ROI-pairs for sim'])
        xlabel('phase lag [rad]')
        ylabel('number of ROI-pairs')
        legend('5 Hz','11 Hz','25 Hz')
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
        export_fig(fullfile(plotDir,[plotFilename '.pdf']), '-transparent', '-r72');

      end
      
      %% plot sim specifics
      if isfield(compareSimExp,'sim')
        if isfield(compareSimExp.sim,'metrics')
          metrics = fieldnames(compareSimExp.sim.spec.metrics);
          for m=1:length(metrics)
            sfigure(1); clf;
            imagesc(compareSimExp.sim.spec.metrics.(metrics{m})(:,:,1,1));
            title(metrics{m})
            colormap hot;
            colorbar
            if length(compareSimExp.sim.spec.dimName)>1
              xlabel(compareSimExp.sim.spec.dimName{2})
              set(gca,'XTick',1:length(compareSimExp.sim.spec.dimLabels{2}));
              set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%10.1f'),compareSimExp.sim.spec.dimLabels{2},'UniformOutput',false))
            end
            ylabel(compareSimExp.sim.spec.dimName{1})
            set(gca,'YTick',1:length(compareSimExp.sim.spec.dimLabels{1}));
            set(gca,'YTickLabel',compareSimExp.sim.spec.dimLabels{1})
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
%             print('-dpng','-r72',fullfile(plotDir,[metrics{m} '.png']))
            export_fig(fullfile(plotDir,[metrics{m} '.pdf']), '-transparent', '-transparent', '-r72');
         end
        end
      end
      
      %% combine pdfs:
      pdfFiles = dir(fullfile(plotDir, '*.pdf'));
      pdfFiles = cellfun(@(x) fullfile(plotDir,x),{pdfFiles.name},'UniformOutput',false);
      if exist([plotDir '.pdf'],'file')
        delete([plotDir '.pdf']);
      end
      append_pdfs([plotDir '.pdf'], pdfFiles{:})
      if p.deletePlotFolder
        rmdir(plotDir,'s');
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
    
    function plotCompareDistWithSC(structure,inSpec,dataSC,caption,filename,plotDir,excludePlotFieldsRegexp)
      if isstruct(structure)
        
        %% load spec if exsits:
        if isfield(structure,'spec')
          inSpec = structure.spec;
        end
        
        %% recursive call:
        fnames = fieldnames(structure);
        for f=1:length(fnames)
          if ~strcmp('spec',fnames{f})
            
            if isempty(regexp(filename,excludePlotFieldsRegexp, 'once'))
              ConnectomeEnvelopeReduce.plotCompareDistWithSC(structure.(fnames{f}), inSpec, dataSC, [caption '.' fnames{f}], [filename '_' fnames{f}], plotDir,excludePlotFieldsRegexp);
            end
            
          end
        end
        
      else
        
        %% plot:
        sfigure(1);
        clf;
        
        if sum(strcmp(inSpec.dimName,'roi'))==2
          % compare with sc_dist:
          trigIds = find(triu(ones(66,66),1));
          if sum(strcmp(inSpec.dimName,'freq'))
            dataSim = squeeze(cell2mat(cellfun(@(x) x(trigIds), num2cell(structure,[1 2]),'UniformOutput',false)));
            dataExp = repmat(dataSC.avg_dist(trigIds),[1 3]);
          else
            dataSim = structure(trigIds);
            dataExp = dataSC.avg_dist(trigIds);
          end
          plot(dataExp,dataSim,'.','MarkerSize', 6);
          xlabel('SC roi distance');
        else
          % compare with sc_roi_size:
          dataSim = structure;
          dataExp = dataSC.avg_roi_size;
          plot(dataExp,dataSim,'.','MarkerSize', 12);
          xlabel('SC roi size');
        end
        
        [RHO,PVAL] = corr(dataExp,dataSim);
        
        ylabel('difference between sim and eeg');
        if isscalar(RHO)
          title([caption ' corr=' num2str(RHO) ' pval=' num2str(PVAL)])
        else
          legend(['5 Hz corr=' num2str(RHO(1,1))],['11 Hz corr=' num2str(RHO(1,2))],['25 Hz corr=' num2str(RHO(1,3))])
          title(caption)
        end
        
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
        
        
%           print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
        export_fig(fullfile(plotDir,[filename '.pdf']), '-transparent', '-r72');
          
      end
      
      
    end
    
    function plotBestParams(structure)
      
      eegSubjDimId = find(strcmp(structure.spec.dimName,'subjEeg'));
      simSubjDimId = find(strcmp(structure.spec.dimName,'subjId'));
      if isempty(eegSubjDimId)
        eegSubjDimId = length(structure.spec.dimName)+1;
        numEegSubj = 1;
      else
        numEegSubj = structure.spec.dimSize(eegSubjDimId);
      end
      if isempty(simSubjDimId)
        simSubjDimId = length(structure.spec.dimName)+1;
        numSimSubj = 1;
      else
        numSimSubj = structure.spec.dimSize(simSubjDimId);
      end
      
      for eegSubjId = 1:numEegSubj
        for simSubjId = 1:numSimSubj
          
          structure.coh
          
        end
      end
      
    end
    
    function dispStructure(structure, caption, fileID, excludeFields)
      % use as follows:
      % fileID = fopen('permTestResults.txt','w');
      % ConnectomeEnvelopeReduce.dispStructure(structure, 'permTest', fileID, {})
      % fclose(fileID);
      
      if isstruct(structure)
        %% recursive call:
        fnames = setdiff(fieldnames(structure),{'spec'});
        for f=1:length(fnames)
          if sum(strcmp(excludeFields,fnames{f}))==0
            ConnectomeEnvelopeReduce.dispStructure(structure.(fnames{f}), [caption '.' fnames{f}], fileID, excludeFields);
          end
        end
      else
        fprintf(fileID,'%s = %12.8f\n',caption,structure);
      end
      
    end
    
    
    function plotStructure(structure, inSpec, caption, filename, plotDir, excludeFields, permutePlotDims, excludePlotFieldsRegexp)
      
      if nargin<7
        permutePlotDims = [];
      end
      
      if nargin<8
        excludePlotFieldsRegexp = '';
      end
      
      if isstruct(structure)
        
        %% load spec if exsits:
        if isfield(structure,'spec')
          inSpec = structure.spec;
        end
        
        %% recursive call:
        fnames = fieldnames(structure);
        for f=1:length(fnames)
          if sum(strcmp(excludeFields,fnames{f}))==0 && ~strcmp('spec',fnames{f})
            
            if isempty(regexp(filename,excludePlotFieldsRegexp, 'once'))
              ConnectomeEnvelopeReduce.plotStructure(structure.(fnames{f}), inSpec, [caption '.' fnames{f}], [filename '_' fnames{f}], plotDir, excludeFields, permutePlotDims, excludePlotFieldsRegexp);
            end
            
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
          if length(siz) == 3
            %% shift third dimension to forth for horizonal plot
            structure = permute(structure,[1 2 4 3]);
            inSpecPerm.dimName{4} = '';
            inSpecPerm.dimLabels{4} = {};
            inSpecPerm.dimName = inSpecPerm.dimName([1 2 4 3]);
            inSpecPerm.dimLabels = inSpecPerm.dimLabels([1 2 4 3]);
          end
          
%           if length(siz) == 3
%             
%             m=round(sqrt(siz(3)));
%             n=ceil(siz(3)/m);
%             
%             sfigure(1); clf;
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
            
            h=sfigure(1);
            clf
            p = panel();
            p.pack(size(structure,3), size(structure,4));
            p.de.margin = 2;
            p.margin = [25 25 20 20];
            p.fontsize = 8;
            p.title(caption);
            for m = 1:size(structure,3)
              for n = 1:size(structure,4)
                p(m, n).select();
                imagesc(structure(:,:,m,n));
                
%                 [nr,nc] = size(structure(:,:,m,n));
%                 pcolor([structure(:,:,m,n) nan(nr,1); nan(1,nc+1)]);
%                 shading flat;
%                 set(gca, 'ydir', 'reverse');
                
                colormap autumn;
                if isnan(minc) || isnan(maxc) || minc==maxc || isinf(maxc) || isinf(minc) 
                  disp('cannot set caxis limits')
                else
                  set(gca,'clim',[minc maxc]);
                end
                
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
                  xticks(xticks>length(inSpecPerm.dimLabels{2})) = [];
                  set(gca,'XTick',xticks);
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
                  yticks(yticks>length(inSpecPerm.dimLabels{1})) = [];
                  set(gca,'YTick',yticks);
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
            saveas(h,fullfile(plotDir,[filename '.fig']));
            
%           else
%             return;
%           end
        elseif siz(2) == 1 && siz(1) == 1
          %% return because only scalar
          return;
        elseif siz(2) == 1
          sfigure(1); clf;
          plot(structure);
          title(caption)
          xlabel(inSpecPerm.dimName{1})
          ylabel(caption)
          set(gca,'XTick',1:size(structure,1));
          set(gca,'XTickLabel',inSpecPerm.dimLabels{1})
          set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
%           print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
          export_fig(fullfile(plotDir,[filename '.pdf']), '-transparent', '-r72');
        else
          h=sfigure(1);clf;
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
%           print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
          saveas(h,fullfile(plotDir,[filename '.fig']));
          export_fig(fullfile(plotDir,[filename '.pdf']), '-transparent', '-r72');
        end
        
      end
      
    end
    
    function plotStructureSubj(structure, inSpec, caption, filename, plotDir, excludeFields)
      
      if isstruct(structure)
        
        %% recursive call:
        fnames = setdiff(fieldnames(structure),{'spec'});
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
          sfigure(1); clf
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
%           print('-dpng','-r72',fullfile(plotDir,[filename '_fix' fixSubj{f} '.png']))
          export_fig(fullfile(plotDir,[filename '_fix' fixSubj{f} '.pdf']), '-transparent', '-r72');
          
        end
      
      end
      
    end
      
    function [tensor, outSpec] = tensorCombineDims(tensor, dimNamesToCombine, method, newDimName, inSpec)
      % method can be 'match', 'nonmatch', 'all', 'trig'
      
      if nargin<5
        inSpec = [];
      end
      
      if isstruct(tensor)
        
        %% load spec if exsits:
        if isfield(tensor,'spec')
          inSpec = tensor.spec;
        end
        
        %% recursive call:
        fnames = setdiff(fieldnames(tensor),{'spec'});
        for f=1:length(fnames)
          [tensor.(fnames{f}), outSpec] = ConnectomeEnvelopeReduce.tensorCombineDims(tensor.(fnames{f}), dimNamesToCombine, method, newDimName, inSpec);
        end
        
        %% update spec from output:
        if isfield(tensor,'spec')
          tensor.spec = outSpec;
        end
        
      else
        
        if iscell(dimNamesToCombine)
          dimId1 = find(strcmp(dimNamesToCombine{1},inSpec.dimName));
          dimId2 = find(strcmp(dimNamesToCombine{2},inSpec.dimName));
          dimId1 = dimId1(1);
          dimId2 = dimId2(end);
        else
          dimIds = find(strcmp(dimNamesToCombine,inSpec.dimName));
          dimId1 = dimIds(1);
          dimId2 = dimIds(2);
        end
        
        
        asCell = num2cell(tensor,[dimId1 dimId2]);
        asCell = cellfun(@(x) squeeze(x), asCell,'UniformOutput', false);
        
        if strcmp(method,'avgMatch')
          
          asCell = cellfun(@(x) mean(diag(x)), asCell,'UniformOutput', false);
          tensor = squeeze(cell2mat(asCell));
          
          outSpec = inSpec;
          outSpec.dimName([dimId1 dimId2]) = [];
          outSpec.dimSize = size(tensor);
          outSpec.dimLabels([dimId1 dimId2]) = [];
          
        else
          if strcmp(method,'match')
            asCell = cellfun(@(x) shiftdim(diag(x),-dimId1+1), asCell,'UniformOutput', false);
          elseif strcmp(method,'nonmatch')
            nonDiagIds=find(~eye([inSpec.dimSize(dimId1) inSpec.dimSize(dimId2)]));
            asCell = cellfun(@(x) shiftdim(x(nonDiagIds),-dimId1+1), asCell,'UniformOutput', false); %#ok<FNDSB>
          elseif strcmp(method,'all')
            asCell = cellfun(@(x) shiftdim(x(:),-dimId1+1), asCell,'UniformOutput', false);
          elseif strcmp(method,'trig')
            trigIds=find(triu(ones([inSpec.dimSize(dimId1) inSpec.dimSize(dimId2)]),1));
            asCell = cellfun(@(x) shiftdim(x(trigIds),-dimId1+1), asCell,'UniformOutput', false); %#ok<FNDSB>
          end
          tensor = permute(cell2mat(asCell),[setdiff(1:length(inSpec.dimSize),dimId2) dimId2]);
          
          outSpec = inSpec;
          outSpec.dimName{dimId1} = newDimName;
          outSpec.dimName(dimId2) = [];
          outSpec.dimSize = size(tensor);
          outSpec.dimLabels{dimId1} = num2cell(1:outSpec.dimSize(dimId1));
          outSpec.dimLabels(dimId2) = [];
        end
        
        
        
        
      end
    end
    
    
    function [structure] = normalizeMat(structure, inSpec, excludeFieldsRegexp )
      
      if isstruct(structure)
        
        %% load spec if exsits:
        if isfield(structure,'spec')
          inSpec = structure.spec;
        end
        
        %% recursive call:
        fnames = setdiff(fieldnames(structure),{'spec'});
        for f=1:length(fnames)
          if isempty(regexp(fnames{f},excludeFieldsRegexp, 'once'))
            structure.(fnames{f}) = ConnectomeEnvelopeReduce.normalizeMat(structure.(fnames{f}), inSpec, excludeFieldsRegexp);
          end
        end
        
      else
        
        eegSubjDim = find(strcmp(inSpec.dimName,'subjEeg'));
        simSubjDim = find(strcmp(inSpec.dimName,'subjId'));
        
        asCell = num2cell(structure,[eegSubjDim simSubjDim]);
        asCell = cellfun(@(x) squeeze(x), asCell,'UniformOutput', false);
        asCell = cellfun(@(x) normMat(x), asCell,'UniformOutput', false);
        
        tmpSize = ones(1,10);
        tmpSize(eegSubjDim) = size(structure,eegSubjDim);
        tmpSize(simSubjDim) = size(structure,simSubjDim);
        
        asCell = cellfun(@(x) reshape(x,tmpSize), asCell,'UniformOutput', false);
        
        structure = cell2mat(asCell);
        
      end
      
    end
    
    
    function [tensorPermTest, outSpec] = evalPermTest(tensor, inSpec, excludeFieldsRegexp )
      
      if nargin<2
        inSpec = [];
      end
      
      if isstruct(tensor)
        
        %% load spec if exsits:
        if isfield(tensor,'spec')
          inSpec = tensor.spec;
        end
        
        %% recursive call:
        fnames = setdiff(fieldnames(tensor),{'spec'});
        for f=1:length(fnames)
          
          if isempty(regexp(fnames{f},excludeFieldsRegexp, 'once'))
            [tensorPermTest.(fnames{f}), outSpec] = ConnectomeEnvelopeReduce.evalPermTest(tensor.(fnames{f}), inSpec, excludeFieldsRegexp );
          end
          
        end
        
        %% update spec from output:
        if isfield(tensor,'spec')
          tensorPermTest.spec = outSpec;
        end
        
      else
        
        eegSubjDim = find(strcmp(inSpec.dimName,'subjEeg'));
        simSubjDim = find(strcmp(inSpec.dimName,'subjId'));
        
        %% do permutation tests for all metrics:
        dimSize = inSpec.dimSize;
        dimSize(end+1:5)=1;
        
        indices = {':',':',':',':',':'};
        iterateOverDims = setdiff(1:5,[eegSubjDim simSubjDim]);
        for iterDim1=1:dimSize(iterateOverDims(1))
          indices{iterateOverDims(1)} = iterDim1;
          for iterDim2=1:dimSize(iterateOverDims(2))
            indices{iterateOverDims(2)} = iterDim2;
            for iterDim3=1:dimSize(iterateOverDims(3))
              indices{iterateOverDims(3)} = iterDim3;
              
              normString = {'orig','norm'};
              for normBefore = 1:2
                RhoMatch = cell(2,1);
                RhoNomat = cell(2,1);
                
                transpString = {'fixEeg','fixSim'};
                for transp=1:2
                  metricVal = squeeze(atanh(tensor(indices{:})));
%                   metricVal = squeeze((tensor(indices{:})));
                  
                  if transp==2
                    metricVal = metricVal';
                  end
                  
                  if normBefore==2
                    metricVal = bsxfun(@rdivide, metricVal, sum(metricVal,1));
                  end
                  
                  %% extract matching and nonmatching matrix:
                  matching = repmat(diag(metricVal),1,size(metricVal,2)-1);
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
                  tensorPermTest.(transpString{transp}).(normString{normBefore}).fractionMatchLarger(iterDim1,iterDim2,iterDim3) = fractionMatchLarger;
                  
                  %% how often is the match also the max in its column
                  
                  
                  %% bootstrap max per column --> p-value ?????
%                   numSubj = length(matching);
%                   numIters = 500;
%                   numSubjPerIter = 2;
%                   for i=1:numIters
%                     selSubj = randi(numSubj,numSubjPerIter);
%                     tmp1 = matching(selSubj);
%                     tmp2 = nonmatch(selSubj,:);
%                     ....
%                   end
                  
                  %% signrank test
                  if sum(isnan(matching(:))) || sum(isnan(nonmatch(:))) || sum(isinf(matching(:))) || sum(isinf(nonmatch(:)))
                    disp(['NaN values in statistical test at iterDim1=' num2str(iterDim1) ' iterDim2=' num2str(iterDim2) ' iterDim3=' num2str(iterDim3) ' '  ])
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).signrank_pval(iterDim1,iterDim2,iterDim3) = NaN;
                  else
                    pval = signrank(matching(:),nonmatch(:));
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).signrank_pval(iterDim1,iterDim2,iterDim3) = pval;
                  end
                  
                  %% do paired t-test with repeated nonmatching variables:
                  RhoMatch{transp} = matching(:);
                  RhoNomat{transp} = nonmatch(:);
                  [~,pval]=ttest(RhoMatch{transp},RhoNomat{transp});
                  tensorPermTest.(transpString{transp}).(normString{normBefore}).repeatMatches(iterDim1,iterDim2,iterDim3) = pval;
                  
                  %% do paired t-test where the nonmatching values are before averaged for each column:
                  nonmatchingMeans = mean(nonmatch,2);
                  [~,pval]=ttest(matching(:,1),nonmatchingMeans);
                  tensorPermTest.(transpString{transp}).(normString{normBefore}).meanNonmatch(iterDim1,iterDim2,iterDim3) = pval;
                  
                  %% do paired t-test where the nonmatching values are before averaged for each column:
                  nonmatchingMedian = median(nonmatch,2);
                  [~,pval]=ttest(matching(:,1),nonmatchingMedian);
                  tensorPermTest.(transpString{transp}).(normString{normBefore}).medianNonmatch(iterDim1,iterDim2,iterDim3) = pval;
                  
                  %% mixed-effects anovan:
                  if ~isreal(metricVal(:)) || sum(isnan(matching(:))) || sum(isnan(nonmatch(:)))
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).anovan_pval_matching(iterDim1,iterDim2,iterDim3) = NaN;
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).anovan_pval_subjDTI(iterDim1,iterDim2,iterDim3) = NaN;
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).anovan_pval_subjEEG(iterDim1,iterDim2,iterDim3) = NaN;
                  else
                    
                    numSubj = size(metricVal,2);
                    groupMatch = eye(numSubj,numSubj);
                    groupSubjDTI = repmat(1:numSubj,[numSubj 1]);
                    groupSubjEEG = repmat(1:numSubj,[numSubj 1])';
                    pval = anovan(metricVal(:), {groupMatch(:) groupSubjDTI(:) groupSubjEEG(:)} , 'random', [2 3], 'varnames', {'matching', 'subjDTI', 'subjEEG'}, 'display','off');
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).anovan_pval_matching(iterDim1,iterDim2,iterDim3) = pval(1);
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).anovan_pval_subjDTI(iterDim1,iterDim2,iterDim3) = pval(2);
                    tensorPermTest.(transpString{transp}).(normString{normBefore}).anovan_pval_subjEEG(iterDim1,iterDim2,iterDim3) = pval(3);
                  end
                  
                end
                
                [~,tensorPermTest.fixBoth(iterDim1,iterDim2,iterDim3)]=ttest([RhoMatch{1}(:); RhoMatch{2}(:)],[RhoNomat{1}(:); RhoNomat{2}(:)]);
                
              end
              
            end
          end
        end
        
        iterateOverDims(iterateOverDims>length(inSpec.dimName)) = [];
        outSpec.dimName = inSpec.dimName(iterateOverDims);
        outSpec.dimSize = size(tensorPermTest.fixSim);
        outSpec.dimLabels = inSpec.dimLabels(iterateOverDims);
        
      end
      
    end
    
    function [tensor, outSpec] = filterTensor(tensor, filterSpec, inSpec )
%       tensor = rand([5,10,2]); or tensor structure
%       filterSpec.subj.Ids = [1:4 7:10];
%       filterSpec.subj.Avg = false;
%       filterSpec.subj.Max = false;
%       filterSpec.day.Ids = [];
%       filterSpec.day.Avg = true;
%       inSpec.dimName = {'day', 'subjEeg', 'cond'};
%       inSpec.dimSize = [5 10 2];
%       inSpec.dimLabels = {'', '', ''};
      
      if nargin<3
        inSpec = [];
      end
      
      if isstruct(tensor)
        
        %% load spec if exsits:
        if isfield(tensor,'spec')
          inSpec = tensor.spec;
        end
        
        %% recursive call:
        fnames = setdiff(fieldnames(tensor),{'spec'});
        for f=1:length(fnames)
          [tensor.(fnames{f}), outSpec] = ConnectomeEnvelopeReduce.filterTensor(tensor.(fnames{f}), filterSpec, inSpec );
        end
        
        %% update spec from output:
        if isfield(tensor,'spec')
          tensor.spec = outSpec;
        end
        
      else
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
          
          if isfield(filterSpec.(filtDimNames{l}),'Avg') && filterSpec.(filtDimNames{l}).Avg
            curDimId = find(strcmp(inSpec.dimName,filtDimNames{l}));
            tensor = mean(tensor,curDimId);
          end
          
          if isfield(filterSpec.(filtDimNames{l}),'Max') &&  filterSpec.(filtDimNames{l}).Max
            curDimId = find(strcmp(inSpec.dimName,filtDimNames{l}));
            tensor = max(tensor,[],curDimId);
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
  
end

