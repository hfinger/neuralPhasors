classdef Syllable < Gridjob
  %GRIDJOBEXAMPLE Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = Syllable(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.Syllable.path = '../Data';
      this.params.Syllable.syllN = 32;
      this.params.Syllable.trainN = 30;
      this.params.Syllable.cvalRuns = 1;
      this.params.Syllable.sampRate = 20000;
      this.params.Syllable.interpolType = 'IIR';
      this.params.Syllable.mfccN = 25;
      this.params.Syllable.invCoeffOrder = 0;
      this.params.Syllable.winsize = 20;
      this.params.Syllable.melFramesN = 64;
      this.params.Syllable.smoothL = 4;
      this.params.Syllable.polyOrder = 3;
      this.params.Syllable.incDer = [1, 1]; %fehlt
      this.params.Syllable.nComp = 10;
      this.params.Syllable.resN = 10;
      this.params.Syllable.specRad = 1.2;
      this.params.Syllable.biasScale = 0.2;
      this.params.Syllable.inpScale = 1.0;
      this.params.Syllable.conn = 1.0;
      this.params.Syllable.gammaPos = 25;
      this.params.Syllable.gammaNeg = 27;
      this.params.Syllable.plotExample = 0;
      this.params.Syllable.targetDir = '';
      this.params.Syllable.scriptsDir = './runSyllClassScripts';
      
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
      source = '/home/student/k/kstandvoss/py/.local/bin/python3.5';
      path = 'cd /net/store/nbp/projects/phasesim/src_kstandvoss/Task1_Recognition/;';
      command = ' runSyllClass.py ';
      options = [this.params.Syllable.path, ' ', num2str(this.params.Syllable.syllN), ' -trainN ', num2str(this.params.Syllable.trainN), ...
          ' -cvalRuns ', num2str(this.params.Syllable.cvalRuns), ' -sampRate ', num2str(this.params.Syllable.sampRate), ' -interpolType ', this.params.Syllable.interpolType, ...
          ' -mfccN ', num2str(this.params.Syllable.mfccN), ' -invCoeffOrder ', num2str(this.params.Syllable.invCoeffOrder), ' -winsize ', num2str(this.params.Syllable.winsize), ...
          ' -melFramesN ', num2str(this.params.Syllable.melFramesN), ' -smoothL ', num2str(this.params.Syllable.smoothL), ' -nComp ', ...
          num2str(this.params.Syllable.nComp), ' -resN ', num2str(this.params.Syllable.resN), ' -specRad ', num2str(this.params.Syllable.specRad), ' -biasScale ', num2str(this.params.Syllable.biasScale), ...
          ' -inpScale ', num2str(this.params.Syllable.inpScale), ' -conn ', num2str(this.params.Syllable.conn), ' -gammaPos ', num2str(this.params.Syllable.gammaPos), ' -gammaNeg ', num2str(this.params.Syllable.gammaNeg), ...
          ' -plotExample ', num2str(this.params.Syllable.plotExample), ' -scriptsDir ', this.params.Syllable.scriptsDir];
      options = [options, ' -targetDir ' , this.params.Syllable.targetDir, '/', num2str(this.currJobid)]; 
      make = ['mkdir ', this.params.Syllable.targetDir, '/', num2str(this.currJobid), ';'];
      call = [path, make, source,  command, options];
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

