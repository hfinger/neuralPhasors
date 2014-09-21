classdef ConnectomeSim < Gridjob
  %ConnectomeSim Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ConnectomeSim(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ConnectomeSim.dataset = 0; % 0=datasimu from Arnaud, 1=SC_Bastian1, 2=dist_and_CI_controls.mat, 3=patients_t1_logCI_mul_20140924_preprocessed, 4=dti_20141209_preprocessed
      this.params.ConnectomeSim.subjId = 1; % or -1 for average subject
      this.params.ConnectomeSim.normRoisizeInterp = []; % 0 = addition, 1 = multiplication
      this.params.ConnectomeSim.dtiDistanceCorrection = false;

      this.params.ConnectomeSim.shuffleSC = false;
      this.params.ConnectomeSim.shuffleD = false;
      this.params.ConnectomeSim.shufflePermutations = [];
      
      this.params.ConnectomeSim.normRowBeforeHomotopic = 0; % 0=no normalization, 1=norm each row, 2=norm the complete matrix,
      this.params.ConnectomeSim.homotopic = 0;
      this.params.ConnectomeSim.roiOutScales = []; % vector of scaling factors for each roi
      this.params.ConnectomeSim.roiOutIds = []; % vector indicating the roiIds to scale
      this.params.ConnectomeSim.normRow = 1; % 0=no normalization, 1=norm each row, 2=norm the complete matrix,
      this.params.ConnectomeSim.model = 'kuramoto'; % 'kuramoto' or 'rate' or 'SAR' or 'WilsonCowan'
      this.params.ConnectomeSim.useNetworkFokkerPlanck = false;

      %params for loading dataset 5:
      this.params.ConnectomeSim.loadSp = 0.6;                          
      this.params.ConnectomeSim.loadHeur = 1;                         
      this.params.ConnectomeSim.loadHscale = 1;  
      
      %params specific for Kuramoto model:
      this.params.ConnectomeSim.approx=false;
      this.params.ConnectomeSim.invertSin=false;
      this.params.ConnectomeSim.f=60;
      this.params.ConnectomeSim.startState = [];
      this.params.ConnectomeSim.saveRelativePhase = false;
     
      %params specific for rate model:
      this.params.ConnectomeSim.tau=20;
      
      %params specific for SAR model:
      this.params.ConnectomeSim.normStd=false;
      
      %params for all models:
      this.params.ConnectomeSim.k=0.8;
      
      %params for all models but SAR model:
      this.params.ConnectomeSim.v=10; % in m/sec
      this.params.ConnectomeSim.delay = 0; % in ms
      this.params.ConnectomeSim.t_max=500;
      this.params.ConnectomeSim.dt=0.0001;
      this.params.ConnectomeSim.sampling=10;
      this.params.ConnectomeSim.sig_n=1.25;
      this.params.ConnectomeSim.d=0;
      this.params.ConnectomeSim.verbose=true;
      
      this.params.ConnectomeSim.statsRemoveInitialT = 0;
            
      this.params.ConnectomeSim.outFilenames = 'results';
      this.params.ConnectomeSim.forceOverwrite = false;

      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});

    end
    
    %% Start: is executed before all individual parameter jobs are started
    function startJob(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%

%       disp('this excecutes before the parallel job is started');
      
      %%%% END EDIT HERE:                                        %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% Run: is executed on the grid for each parameter individually (in parallel)
    % this function is called from Gridjob-class and executes the parallel job
    function run(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: implement the algorithm here %%%%
      
      param = this.params.ConnectomeSim;
      
      paths = dataPaths( );
      
      savepath = fullfile(this.workpath,this.params.ConnectomeSim.outFilenames);
      if param.forceOverwrite || ~exist(fullfile(savepath,[num2str(this.currJobid) 'FC.mat']),'file')
        
        if param.dataset==0
          load(fullfile(paths.databases,'datasimu.mat'));
        elseif param.dataset==1
          
          load(fullfile(paths.databases,'SC_Bastian','Stuct_connectivity_for_Holger.mat'));
          
          if param.subjId==-1
            SC = mean(cell2mat(reshape(SC,[1 1 numel(SC)])),3);
            D = mean(cell2mat(reshape(D,[1 1 numel(D)])),3);
          else
            SC = SC{param.subjId};
            D = D{param.subjId};
          end
          
          SC(isnan(SC)) = 0;
          SC = SC + SC';
          
          D(D==0)=Inf;
          D(isinf(D))=min(D(:));
          D = 1./D;
          D(isnan(D)) = 0;
          D = D + D';
          D = D * 230; % approx to scale to mm
        elseif param.dataset==2 || param.dataset==4
          if param.dataset==2
            ci = load(fullfile(paths.databases,'SC_Bastian','dist_and_CI_controls_preprocessed.mat'));
          elseif param.dataset==4
            ci = load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));
          end
          
          if ~isempty(param.normRoisizeInterp)
            if length(param.subjId)>1
              avg_ci = zeros(size(ci.samples{1}));
              counter = 0;
              for ss = 1:length(param.subjId)
                if ~isempty(ci.samples{param.subjId(ss)})
                  samples = ci.samples{param.subjId(ss)};
                  roisize = ci.roisize{param.subjId(ss)};
                  roisizeMul = roisize * roisize';
                  roisizeAdd = bsxfun(@plus, roisize, roisize');
                  avg_ci = avg_ci + samples ./ (param.normRoisizeInterp * roisizeMul + (1-param.normRoisizeInterp) * roisizeAdd);
                  counter = counter + 1;
                end
              end
              SC = avg_ci / counter;
              D = ci.avg_dist;
            elseif param.subjId==-1
              avg_ci = zeros(size(ci.samples{1}));
              counter = 0;
              for subjId = 1:length(ci.samples)
                if ~isempty(ci.samples{subjId})
                  samples = ci.samples{subjId};
                  roisize = ci.roisize{subjId};
                  roisizeMul = roisize * roisize';
                  roisizeAdd = bsxfun(@plus, roisize, roisize');
                  avg_ci = avg_ci + samples ./ (param.normRoisizeInterp * roisizeMul + (1-param.normRoisizeInterp) * roisizeAdd);
                  counter = counter + 1;
                end
              end
              SC = avg_ci / counter;
              D = ci.avg_dist;
            else
              samples = ci.samples{param.subjId};
              roisize = ci.roisize{param.subjId};
              roisizeMul = roisize * roisize';
              roisizeAdd = bsxfun(@plus, roisize, roisize');
              SC = samples ./ (param.normRoisizeInterp * roisizeMul + (1-param.normRoisizeInterp) * roisizeAdd);
              D = ci.dist{param.subjId};
            end
          else
            if param.subjId==-1
              SC = ci.avg_ci;
              D = ci.avg_dist;
            else
              if length(param.subjId)>1
                SC = nanmean(cell2mat(permute(ci.ci(param.subjId),[1 3 2])),3);
                D = nanmean(cell2mat(permute(ci.dist(param.subjId),[1 3 2])),3);
              else
                SC = ci.ci{param.subjId};
                D = ci.dist{param.subjId};
              end
            end
          end
          
          SC(logical(tril(ones(size(SC))))) = 0;
          SC(isnan(SC)) = 0;
          SC = SC + SC';
          
          D(logical(tril(ones(size(D))))) = 0;
          D(isnan(D)) = 0;
          D = D + D'; % scale unclear should be scaled to mm
          
        elseif param.dataset==3
          ci = load(fullfile(paths.databases,'SC_Bastian','patients_t1_logCI_mul_20140924_preprocessed.mat'));
          
          if ~isempty(param.normRoisizeInterp)
            samples = ci.samples{param.subjId};
            roisize = ci.roisize{param.subjId};
            roisizeMul = roisize * roisize';
            roisizeAdd = bsxfun(@plus, roisize, roisize');
            SC = samples ./ (param.normRoisizeInterp * roisizeMul + (1-param.normRoisizeInterp) * roisizeAdd);
            D = ci.dist{param.subjId};
          else
            SC = ci.ci{param.subjId};
            D = ci.dist{param.subjId};
          end
          
          SC(isnan(SC)) = 0;
          %         SC = SC + SC';
          
          D(isnan(D)) = 0;
          %         D = D + D'; % scale unclear should be scaled to mm
        
        elseif param.dataset==5
          
          jb = load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics','temp_ConnectomeMetrics','jobDesc.mat'));
          % xy.paramComb % column-wise representation of jobIDs 
          % xy.variableParams{variable,1}{1,2};
          
          getId = find(sum(ismember(cell2mat(jb.paramComb),[param.loadSp; param.loadHeur; param.loadHscale]),1) == length(jb.variableParams),1);    

          ci = load(fullfile(paths.workdir,'pebel','20150414_SAR_Metrics','results',strcat(num2str(getId),'SC.mat')));        
          SC = ci.hSC;
          D  = ci.hMetr.perConn.euclDist; % distance matrix for distance correction, see param ConnectomeSim.dtiDistanceCorrection
        end
        
        if param.shuffleSC || param.shuffleD
          trigIds = find(triu(ones(size(SC)),1));
          
          if ~isempty(param.shufflePermutations)
            p = param.shufflePermutations;
          else
            p = randperm(length(trigIds(:)),length(trigIds(:)));
          end
          
          if param.shuffleSC
            tmp = SC(trigIds);
            tmp= tmp(p);
            SC = zeros(size(SC));
            SC(trigIds) = tmp;
            SC = SC + SC';
          end
          if param.shuffleD
            tmp = D(trigIds);
            tmp= tmp(p);
            D = zeros(size(D));
            D(trigIds) = tmp;
            D = D + D';
          end
        end
        
        
        if param.dtiDistanceCorrection
          SC = SC.*D;
        end
        
        if param.normRowBeforeHomotopic==1
          SC = bsxfun(@rdivide,SC,sum(SC,2));
        elseif param.normRowBeforeHomotopic==2
          SC = size(SC,1) * SC / sum(SC(:));
        end
        
        if param.homotopic
          SC = SC + param.homotopic * diag(ones(size(SC,1)/2,1),size(SC,1)/2) + param.homotopic * diag(ones(size(SC,1)/2,1),-size(SC,1)/2);
        end
        
        numRois = size(SC,1);
        if ~isempty(param.roiOutScales)
          
          if ~isempty(param.roiOutIds)
            roiOutScales = ones(1,numRois);
            roiOutScales(param.roiOutIds) = param.roiOutScales;
          else
            roiOutScales = param.roiOutScales;
          end
          
          roiOutScales = repmat(reshape(roiOutScales,[1 numRois]),[numRois 1]);
          SC = SC .* roiOutScales;
          
        end
        
        if param.normRow==1
          C = bsxfun(@rdivide,SC,sum(SC,2));
        elseif param.normRow==2
          C = size(SC,1) * SC / sum(SC(:));
        else
          C = SC;
        end
        
        if param.delay
          %% incorporate the fixed time delay offset into the distance matrix
          % 1. convert distance (in mm) to delay (in ms) via velocity:
          timeDelays = D / param.v;
          % 2. add constant delay offset (in ms)
          timeDelays = timeDelays + param.delay;
          % 3. convert back to distance (in mm):
          D = timeDelays * param.v;
        end
        
        savepath = fullfile(this.workpath,this.params.ConnectomeSim.outFilenames);
        if ~exist(savepath,'dir')
          mkdir(savepath);
        end
        savefilename = fullfile(savepath,num2str(this.currJobid));
        
        if strcmp(param.model,'kuramoto')
          if param.useNetworkFokkerPlanck
            phase = Network_FokkerPlanckKuramoto(C,D,param.k,param.f,param.v,param.t_max,param.dt,param.sampling,param.sig_n,param.d,param.verbose,param.approx,param.invertSin,param.startState);
          else
            phase = runKuramoto(C,D,param.k,param.f,param.v,param.t_max,param.dt,param.sampling,param.sig_n,param.d,param.verbose,param.approx,param.invertSin,param.startState);
          end
          
          cutAt = max(1,round(param.statsRemoveInitialT/(param.dt*param.sampling)));
          OrderParam = abs(sum(exp(1i*phase(:,cutAt:end)),1)/size(phase,1));
          meanOrderParam = mean(OrderParam);
          stdOrderParam = std(OrderParam);
          save([savefilename 'KuramotoStats.mat'],'OrderParam','meanOrderParam','stdOrderParam');
          save([savefilename 'SimResult.mat'],'phase','param');
          
          FCsimNoBold = corr(sin(phase)');
          
          if param.saveRelativePhase
            phase = bsxfun(@minus, phase, 2 * pi * param.f * param.dt * param.sampling * (1:size(phase,2)));
            save([savefilename 'SimResultRelativePhase.mat'],'phase','param');
          end
          
        elseif strcmp(param.model,'WilsonCowan')
          
          [u, v] = runWilsonCowan(C,D,param.k,param.f,param.v,param.t_max,param.dt,param.sampling,param.sig_n,param.d,param.verbose,param.approx,param.invertSin,param.startState);
          
          
        elseif strcmp(param.model,'Breakspear')
          
          
          
        elseif strcmp(param.model,'rate')
          if param.useNetworkFokkerPlanck
            rate = Network_FokkerPlanck(C,D,param.k,param.tau,param.v,param.t_max,param.dt,param.sampling,param.sig_n,param.d,param.verbose);
          else
            rate = runRatemodel(C,D,param.k,param.tau,param.v,param.t_max,param.dt,param.sampling,param.sig_n,param.d,param.verbose);
          end
          save([savefilename 'SimResult.mat'],'rate','param');
          
          FCsimNoBold = corr(rate');
          
        elseif strcmp(param.model,'SAR')
          
          tmp = eye(size(C))-param.k*C;
          FCsimNoBold = inv(tmp) * inv(tmp');
          
          if this.params.ConnectomeSim.normStd
            stdev=sqrt(diag(FCsimNoBold));
            FCsimNoBold = FCsimNoBold./(stdev*stdev');
          end
          
        end
        
        FCsimNoBoldVals = FCsimNoBold(find(tril(ones(size(FCsimNoBold)),-1)));
        if exist('FC','var')
          FCVals = FC(find(tril(ones(size(FC)),-1)));
          [corrFCvsFCsimNoBold,pvalFCvsFCsimNoBold] = corr(FCsimNoBoldVals(:),FCVals(:));
          save([savefilename 'FC.mat'],'FCsimNoBold','corrFCvsFCsimNoBold','pvalFCvsFCsimNoBold');
        else
          save([savefilename 'FC.mat'],'FCsimNoBold');
        end
        
      end
      
      if isfield(param,'ConnectomeEnvelope')
        for ceId=1:length(param.ConnectomeEnvelope)
          ceParams.ConnectomeEnvelope = param.ConnectomeEnvelope(ceId);
          ceParams.ConnectomeEnvelope.inFileRates = fullfile(this.params.ConnectomeSim.outFilenames,[num2str(this.currJobid) 'SimResult.mat']);
          ce = ConnectomeEnvelope(ceParams);
          ce.workpath = this.workpath;
          ce.currJobid = this.currJobid;
          ce.run();
        end
      end
      
      
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
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function plotJob(this)
      
      %% plot instantanous frequency of all neurons:
      figure;
      imagesc(diff(Y2/(2*pi),[],2))
      
      plot(OrderParam)
      
      
      
    end
    
  end
  
end

