classdef ComplexAutoencoder < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ComplexAutoencoder(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ComplexAutoencoder.hiddenN = 50;
      this.params.ComplexAutoencoder.kernelSize = 7;
      this.params.ComplexAutoencoder.trainingSteps = 5000;
      this.params.ComplexAutoencoder.batchsize = 64;
      this.params.ComplexAutoencoder.learningRate = 0.05;
      this.params.ComplexAutoencoder.weightDecay = 5E-5;
      this.params.ComplexAutoencoder.learningRateDecay = 5E-7;
      this.params.ComplexAutoencoder.noiseRate = 0.1;
      this.params.ComplexAutoencoder.momentum = 0;
      this.params.ComplexAutoencoder.coefL1 = 1e-5;
      this.params.ComplexAutoencoder.gpu = 0;   
%      this.comments.GridjobExample.param1 = 'this is some test number parameter';
%      this.comments.GridjobExample.param2 = 'this is some test boolean parameter';
%      this.comments.GridjobExample.param3 = 'this is some test string parameter';
      
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
      
      disp('Start Jobs')
      export = 'export PATH=:$PATH:/net/store/nbp/projects/phasesim/local/bin ;';
      path = 'th /net/store/nbp/projects/phasesim/src_kstandvoss/lua/autoencoder/ComplexAutoencoder/src/main.lua ';
      options = [' -f ', num2str(this.params.ComplexAutoencoder.trainingSteps), ' -r ', num2str(this.params.ComplexAutoencoder.learningRate) ...
          , ' -b ', num2str(this.params.ComplexAutoencoder.batchsize), ' -d ', num2str(this.params.ComplexAutoencoder.weightDecay) ...
          , ' -h ', num2str(this.params.ComplexAutoencoder.hiddenN), ' -k ', num2str(this.params.ComplexAutoencoder.kernelSize) ...
          , ' -n ', num2str(this.params.ComplexAutoencoder.noiseRate), ' -l ', num2str(this.params.ComplexAutoencoder.coefL1)...
          , ' -m ', num2str(this.params.ComplexAutoencoder.momentum)];
      options = [options, ' - w ' , this.params.Gridjob.workpath, ' -i ' this.currJobid]; 
      if this.params.ComplexAutoencoder.gpu 
          options = [options, ' -u ']; 
      end
      call = [export, path, options];
      system(call);

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

