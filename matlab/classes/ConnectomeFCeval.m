classdef ConnectomeFCeval < Gridjob
  %ConnectomeFCeval Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ConnectomeFCeval(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ConnectomeFCeval.inFileRates = 'results';
      this.params.ConnectomeFCeval.t_rm = 20;
      this.params.ConnectomeFCeval.outFilenames = 'boldSignal';
      
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
      
      input = load(fullfile(this.workpath, this.params.ConnectomeFCeval.inFileRates));
      if strcmp(input.param.model,'kuramoto')
        if isfield(input,'phase')
          phase = input.phase;
          rate = sin(phase);
        else
          rate = input.Y;
        end
      else
        rate = input.rate;
      end
      
      if isfield(input,'param')
        dt = input.param.dt;
        sampling = input.param.sampling;
      else
        dt = 0.0001;
        sampling = 10;
      end
      
      % COMPUTE BOLD SIGNAL ==============================
      % using the nonlinear balloon-windkessel model...
      disp('beginning bold calculation ...');
      tic;
      Ybold = mBOLDs(rate,dt*sampling,this.params.ConnectomeFCeval.t_rm);
      time_bold = toc;
      
      % COMPUTE SOME BASIC BOLD SIGNAL ANALYSES ==========
      % settings for bold averages
      xsec = 2000; xgap = 500;        % in msec, window size and spacing
      
      T = size(Ybold,2)*dt*sampling*1000;
      t0 = 1:xgap:T-xsec+1; t0 = ceil(t0/(dt*sampling*1000));
      te = xsec:xgap:T; te = ceil(te/(dt*sampling*1000));
      
      % compute bold averages
      R = size(Ybold,1);
      Ybold_w = zeros(R,length(t0));
      for w=1:length(t0), Ybold_w(:,w) = mean(Ybold(:,t0(w):te(w)),2); end
      
      % remove NaNs, get average bold signal over whole brain, and regress out
      Ybold_w(isnan(Ybold_w)) = 0;
      Ybold_w_mean = mean(Ybold_w,1);
      Ybold_w_reg = zeros(R,length(t0));
      if ~isempty(Ybold_w), 
        for i=1:R, 
          [~,~,Ybold_w_reg(i,:)] = regress(Ybold_w(i,:)',Ybold_w_mean'); 
        end, 
      end
      
      % calc FC:
      FCsim = corr(Ybold_w_reg');
      
      % compare FC simulated and measured:
      paths = dataPaths( );
      load(fullfile(paths.databases,'datasimu.mat'));
      
      FCsimVals = FCsim(find(tril(ones(size(FCsim)),-1)));
      FCVals = FC(find(tril(ones(size(FC)),-1)));
      [corrFCvsFCsim,pvalFCvsFCsim] = corr(FCsimVals(:),FCVals(:));
      
      savepath = fullfile(this.workpath,this.params.ConnectomeFCeval.outFilenames);
      if this.numJobs > 1
        savepath = fullfile(savepath,num2str(this.currJobid));
      end
      mkdir(savepath);
      save(savepath,'Ybold_w_reg','Ybold_w','FCsim','corrFCvsFCsim','pvalFCvsFCsim');
      
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

