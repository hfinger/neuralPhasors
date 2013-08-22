classdef GridToyKuramoto < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = GridToyKuramoto(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.GridToyKuramoto.sim.k=12;
      this.params.GridToyKuramoto.sim.f=40;
      this.params.GridToyKuramoto.sim.v=10;
      this.params.GridToyKuramoto.sim.t_max=10;
      this.params.GridToyKuramoto.sim.dt=0.0001;
      this.params.GridToyKuramoto.sim.sampling=10;
      this.params.GridToyKuramoto.sim.sig_n=0;
      this.params.GridToyKuramoto.sim.d=0;
      this.params.GridToyKuramoto.sim.verbose=true;
      this.params.GridToyKuramoto.sim.approx=false;
      this.params.GridToyKuramoto.sim.invertSin=false;
      
      this.params.GridToyKuramoto.network.weightedNetwork=false;
      this.params.GridToyKuramoto.network.binaryNetwork=true;
      this.params.GridToyKuramoto.network.N=5;
      this.params.GridToyKuramoto.network.k_intraCluster = 0.9;
      this.params.GridToyKuramoto.network.delay_intraClust = 2; %in ms
      this.params.GridToyKuramoto.network.numCluster = 3;
      this.params.GridToyKuramoto.network.k_interCluster = 0.1;
      this.params.GridToyKuramoto.network.delay_interClust = 5; %in ms
      
      this.params.GridToyKuramoto.env.t_rm=4;
      
      this.params.GridToyKuramoto.env.sigBandpass(1).Fst1 = 5.5;
      this.params.GridToyKuramoto.env.sigBandpass(1).Fp1 = 6;
      this.params.GridToyKuramoto.env.sigBandpass(1).Fp2 = 22;
      this.params.GridToyKuramoto.env.sigBandpass(1).Fst2 = 22.5;
      this.params.GridToyKuramoto.env.sigBandpass(1).Ast1 = 40;
      this.params.GridToyKuramoto.env.sigBandpass(1).Ap = 1;
      this.params.GridToyKuramoto.env.sigBandpass(1).Ast2 = 40;
      
      this.params.GridToyKuramoto.env.sigBandpass(2).Fst1 = 29.5;
      this.params.GridToyKuramoto.env.sigBandpass(2).Fp1 = 30;
      this.params.GridToyKuramoto.env.sigBandpass(2).Fp2 = 48;
      this.params.GridToyKuramoto.env.sigBandpass(2).Fst2 = 48.5;
      this.params.GridToyKuramoto.env.sigBandpass(2).Ast1 = 40;
      this.params.GridToyKuramoto.env.sigBandpass(2).Ap = 1;
      this.params.GridToyKuramoto.env.sigBandpass(2).Ast2 = 40;
      
      this.params.GridToyKuramoto.env.envLowpass.Fp = 0.5;
      this.params.GridToyKuramoto.env.envLowpass.Fst = 1;
      this.params.GridToyKuramoto.env.envLowpass.Ap = 1;
      this.params.GridToyKuramoto.env.envLowpass.Ast = 40;
      
      
      
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
      
      disp('some job parameters:')
      disp(this.workpath);
      disp(this.temppath);
      disp(this.resultpath);
      disp(this.currJobid);
      
      [ C, D ] = genNetwork( this.params.GridToyKuramoto.network );
      [ simResult ] = runKuramotoSim( this.params.GridToyKuramoto.sim, C , D );
      [ simEval ] = calcEnvFC( this.params.GridToyKuramoto.env, simResult );
      save(fullfile(this.workpath,['sim' num2str(this.currJobid)]));
      
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

