classdef ProbtrackX < Gridjob
  %ProbtrackX Summary of this class goes here
  %   Detailed explanation goes here
  
  properties

  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ProbtrackX(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ProbtrackX.split = cell2mat(1:53);
      this.params.ProbtrackX.numberPerSplit = 1000;
      this.params.ProbtrackX.numberRegions = 52442;
      
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
      disp(this.params.ProbtrackX.param1);
      
      compl_fs_mask = load_untouch_nii('/net/store/nbp/phasesim/workdir/Arushi/20150323 - single voxel run with term mask/compl_fs_mask.nii');
      new_tract_space = load('/net/store/nbp/phasesim/workdir/Arushi/20150323 - single voxel run with term mask/new_tract_space.mat');
      first_voxel = ((this.params.ProbtrackX.split-1) * 1000) +1;
      
      if this.params.ProbtrackX.split ~= 53
        for i=first_voxel:first_voxel+this.params.ProbtrackX.numberPerSplit-1
                x = new_tract_space(i,1:3);
  
                if compl_fs_mask.img(x(1),x(2),x(3)) == 1
    
                    compl_fs_mask.img(x(1),x(2),x(3)) = 0;
                    save_untouch_nii(compl_fs_mask, ['/net/store/nbp/phasesim/workdir/Arushi/20150323 - single voxel run with term mask/voxel' num2str(i)]);
                    compl_fs_mask.img(x(1),x(2),x(3)) = 1;
                    
                    %include call for probtrackx2
                   
                else
                    disp('voxel', i, 'is not populated in the fs mask');
                end
        end
      else
        for i = first_voxel:size(new_tract_space)
                x = new_tract_space(i,1:3);
  
                if compl_fs_mask.img(x(1),x(2),x(3)) == 1 
    
                    compl_fs_mask.img(x(1),x(2),x(3)) = 0;
                    save_untouch_nii(compl_fs_mask, ['/net/store/nbp/phasesim/workdir/Arushi/20150323 - single voxel run with term mask/voxel' num2str(i)]);
                    compl_fs_mask.img(x(1),x(2),x(3)) = 1;
                    
                    %include call for probtrackx2
                    
                else
                    disp('voxel', i, 'is not populated in the fs mask');
                end
        end
      end
      
      
         
      save(fullfile(this.workpath,[this.currJobid '.mat']));
      
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

