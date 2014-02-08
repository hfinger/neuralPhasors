classdef ConnectomeEnvelope < Gridjob
  %ConnectomeEnvelope Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    Fs
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ConnectomeEnvelope(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ConnectomeEnvelope.inFileRates = 'results';
      this.params.ConnectomeEnvelope.downsamplingFactor = 1;
      this.params.ConnectomeEnvelope.source_t_start = 19;
      this.params.ConnectomeEnvelope.source_t_end = Inf;
      this.params.ConnectomeEnvelope.saveSamples_t_start = -Inf;
      this.params.ConnectomeEnvelope.saveSamples_t_end = Inf;
      this.params.ConnectomeEnvelope.env_t_start = 1;
      this.params.ConnectomeEnvelope.env_t_end = Inf;
      this.params.ConnectomeEnvelope.applyLeadField = false;
      this.params.ConnectomeEnvelope.applyLeadFieldMethod = 'sum'; %or 'eucl'
      this.params.ConnectomeEnvelope.filtermethod = 'butter'; %or equiripple
      
      if length(this.params.ConnectomeEnvelope.sigBandpass)>1
        for b=1:length(this.params.ConnectomeEnvelope.sigBandpass)
          this.params.ConnectomeEnvelope.sigBandpass(b).Fst1 = 10;
          this.params.ConnectomeEnvelope.sigBandpass(b).Fp1 = 10.5;
          this.params.ConnectomeEnvelope.sigBandpass(b).Fp2 = 21.5;
          this.params.ConnectomeEnvelope.sigBandpass(b).Fst2 = 22;
          this.params.ConnectomeEnvelope.sigBandpass(b).Ast1 = 40;
          this.params.ConnectomeEnvelope.sigBandpass(b).Ap = 1;
          this.params.ConnectomeEnvelope.sigBandpass(b).Ast2 = 40;
        end
      else
        this.params.ConnectomeEnvelope.sigBandpass.Fst1 = 10;
        this.params.ConnectomeEnvelope.sigBandpass.Fp1 = 10.5;
        this.params.ConnectomeEnvelope.sigBandpass.Fp2 = 21.5;
        this.params.ConnectomeEnvelope.sigBandpass.Fst2 = 22;
        this.params.ConnectomeEnvelope.sigBandpass.Ast1 = 40;
        this.params.ConnectomeEnvelope.sigBandpass.Ap = 1;
        this.params.ConnectomeEnvelope.sigBandpass.Ast2 = 40;
      end
      
      this.params.ConnectomeEnvelope.envLowpass.Fp = 0.5;
      this.params.ConnectomeEnvelope.envLowpass.Fst = 1;
      this.params.ConnectomeEnvelope.envLowpass.Ap = 1;
      this.params.ConnectomeEnvelope.envLowpass.Ast = 40;
      
      this.params.ConnectomeEnvelope.outFilenames = 'connEnv';
      this.params.ConnectomeEnvelope.saveSourceRate = false;
      this.params.ConnectomeEnvelope.saveSourcePhase = false;
      this.params.ConnectomeEnvelope.saveSourceBP = false;
      this.params.ConnectomeEnvelope.saveEnvSig = false;
      this.params.ConnectomeEnvelope.saveEnvLP = false;
      
      this.params.ConnectomeEnvelope.deleteSimResult = false;
      
      
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
      
      p = this.params.ConnectomeEnvelope;
      
      input = load(fullfile(this.workpath, this.params.ConnectomeEnvelope.inFileRates));
      if isfield(input,'param')
        this.Fs = 1 / (input.param.sampling * input.param.dt);
      else
        this.Fs = 1000;
      end
      
      savepath = fullfile(this.workpath,this.params.ConnectomeEnvelope.outFilenames);
      if ~exist(savepath,'dir')
        mkdir(savepath);
      end

      if strcmp(input.param.model,'kuramoto')
        if isfield(input,'phase')
          phase = this.extracTimeWindows(input.phase,p.source_t_start,p.source_t_end);
          phase = phase{1};
          rate = sin(phase);
        else
          rate = this.extracTimeWindows(input.Y,p.source_t_start,p.source_t_end);
          rate = rate{1};
        end
      else
        rate = this.extracTimeWindows(input.rate,p.source_t_start,p.source_t_end);
        rate = rate{1};
      end
      
      if p.applyLeadField
        paths = dataPaths();
        lf = load(fullfile(paths.databases,'SC_Bastian','20141006_lf_all.mat'));
        lf = lf.lf_all{input.param.subjId};
        lf = cell2mat(permute(lf,[1 3 2]));
        
        % first just use sum of 3 dimensions:
        if strcmp(p.applyLeadFieldMethod,'sum')
          lf = sum(lf,2);
          lf = permute(lf,[1 3 2]);
          rate = lf*rate;
        elseif strcmp(p.applyLeadFieldMethod,'eucl')
          lf = sqrt(sum(lf.^2,2));
          lf = permute(lf,[1 3 2]);
          rate = lf*rate;
          rate = bsxfun(@minus, rate, mean(rate,2));
        elseif strcmp(p.applyLeadFieldMethod,'maxpower')
          % rate = [66 x 50000]
          % lf = [63 x 3 x 66]
          for o=1:3
            tmp = squeeze(lf(:,o,:));
            rate2(:,:,o) = tmp * rate;
          end
          sigPower = squeeze(rms(rate2,2));
          [C,I] = max(sigPower,[],2);
          
          for chan=1:size(rate2,1)
            rate(chan,:) = rate2(chan,:,I(chan));
          end
          
        end
        
      end
      
      if p.downsamplingFactor~=1
        rate = downsample(rate', p.downsamplingFactor)';
        phase = downsample(phase', p.downsamplingFactor)';
        this.Fs = this.Fs / p.downsamplingFactor;
      end
      
      if p.saveSourceRate
        Y = this.extracTimeWindows(rate,p.saveSamples_t_start,p.saveSamples_t_end); %#ok<NASGU>
        save([savepath '/sourceRate_job' num2str(this.currJobid)],'Y');
      end
      
      if p.saveSourcePhase
        Y = this.extracTimeWindows(phase,p.saveSamples_t_start,p.saveSamples_t_end); %#ok<NASGU>
        save([savepath '/sourcePhase_job' num2str(this.currJobid)],'Y');
      end

      FCsim = cell(1,length(p.sigBandpass));
      for b=1:length(p.sigBandpass)
        
        disp('beginning band-pass filter ...');
        bp = p.sigBandpass(b);
        dBandPass = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',bp.Fst1,bp.Fp1,bp.Fp2,bp.Fst2,bp.Ast1,bp.Ap,bp.Ast2,this.Fs);
        HdBandPass = design(dBandPass,'butter');
        sourceBP = zeros(size(rate));
        for k=1:size(rate,1)
          sourceBP(k,:) = filter(HdBandPass,rate(k,:));
        end
        if p.saveSourceBP
          Y = this.extracTimeWindows(sourceBP,p.saveSamples_t_start,p.saveSamples_t_end);
          save([savepath '/sourceBP_band' num2str(b) '_job' num2str(this.currJobid)],'Y');
        end
        
        disp('calc analytic signal ...');
        analyticSig = zeros(size(rate));
        for k=1:size(rate,1)
          analyticSig(k,:) = hilbert(sourceBP(k,:));
        end
        
        disp('calc PLV ...');
        phaseSig = analyticSig./abs(analyticSig);
        phaseSigForPLV = this.extracTimeWindows(phaseSig,p.env_t_start,p.env_t_end);
        PLVsim{b}.bp = bp;
        for s=1:length(phaseSigForPLV)
          PLVsim{b}.PLV{s} = zeros(size(phaseSig,1),size(phaseSig,1));
          for k=1:size(phaseSig,1)
            PLVsim{b}.PLV{s}(k,:) = abs(mean( bsxfun(@times, conj(phaseSigForPLV{s}(k,:)), phaseSigForPLV{s}) , 2 ));
          end
        end
        
        % phaseSig = [66 x 50000] = [ ROI x TIME ] 
        
        disp('calc imag coherence ...');
        analyticSigForICOH = this.extracTimeWindows(analyticSig,p.env_t_start,p.env_t_end);
        ICOHsim{b}.bp = bp;
        for s=1:length(analyticSigForICOH)
          ICOHsim{b}.autospectrum{s} = mean( conj(analyticSigForICOH{s}) .* analyticSigForICOH{s} , 2 );
          ICOHsim{b}.crossspectrum{s} = zeros(size(analyticSigForICOH{s},1),size(analyticSigForICOH{s},1));
          ICOHsim{b}.meanImagCrossSpec{s} = zeros(size(analyticSigForICOH{s},1),size(analyticSigForICOH{s},1));
          ICOHsim{b}.meanAbsImagCrossSpec{s} = zeros(size(analyticSigForICOH{s},1),size(analyticSigForICOH{s},1));
          ICOHsim{b}.pli{s} = zeros(size(analyticSigForICOH{s},1),size(analyticSigForICOH{s},1));
          for k=1:size(analyticSigForICOH{s},1)
            X = bsxfun(@times, conj(analyticSigForICOH{s}(k,:)), analyticSigForICOH{s});
            ICOHsim{b}.crossspectrum{s}(k,:) = mean( X , 2 );
            ICOHsim{b}.meanImagCrossSpec{s}(k,:) = mean(imag(X),2);
            ICOHsim{b}.meanAbsImagCrossSpec{s}(k,:) = mean(abs(imag(X)),2);
            ICOHsim{b}.pli{s}(k,:) = mean(sign(imag(X)),2);
          end
          ICOHsim{b}.wpli{s} = ICOHsim{b}.meanImagCrossSpec{s}./ICOHsim{b}.meanAbsImagCrossSpec{s};          
          ICOHsim{b}.coherency{s} = ICOHsim{b}.crossspectrum{s} ./ sqrt( ICOHsim{b}.autospectrum{s} * ICOHsim{b}.autospectrum{s}' );
          ICOHsim{b}.coherence{s} = abs(ICOHsim{b}.coherency{s});
          ICOHsim{b}.icoh{s} = imag(ICOHsim{b}.coherency{s});
        end
        
        disp('calc envelope signal ...');
        envSig = abs(analyticSig);
        if p.saveEnvSig
          Y = this.extracTimeWindows(envSig,p.saveSamples_t_start,p.saveSamples_t_end);
          save([savepath '/envSig_band' num2str(b) '_job' num2str(this.currJobid)],'Y');
        end
        
        disp('beginning low-pass filter ...');
        lp = p.envLowpass;
        dLowPass = fdesign.lowpass('Fp,Fst,Ap,Ast',lp.Fp,lp.Fst,lp.Ap,lp.Ast,this.Fs);
        HdLowPass = design(dLowPass,'butter');
        envLP = zeros(size(envSig));
        for k=1:size(envSig,1)
          envLP(k,:) = filter(HdLowPass,envSig(k,:));
        end
        if p.saveEnvLP
          Y = this.extracTimeWindows(envLP,p.saveSamples_t_start,p.saveSamples_t_end);
          save([savepath '/envLP_band' num2str(b) '_job' num2str(this.currJobid)],'Y');
        end
        
        envLPforFC = this.extracTimeWindows(envLP,p.env_t_start,p.env_t_end);
        FCsim{b}.bp = bp;
        for s=1:length(envLPforFC)
          FCsim{b}.FC{s} = corr(envLPforFC{s}');
        end
      end
      
      save([savepath '/FCsimJob' num2str(this.currJobid)],'FCsim','PLVsim','ICOHsim');
      
      if p.deleteSimResult
        delete(fullfile(this.workpath, this.params.ConnectomeEnvelope.inFileRates))
      end
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    function savevar = extracTimeWindows(this,fullvar,t_start,t_end)
      savevar = cell(1,length(t_start));
      for s=1:length(t_start)
        cutStart = max(1,round(t_start(s)*this.Fs));
        cutEnd = min(size(fullvar,2),round(t_end(s)*this.Fs));
        savevar{s} = fullvar(:,cutStart:cutEnd);
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
  
end

