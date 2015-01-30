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
      
      this.params.ProbtrackX.split = num2cell(1:53);
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
      
%       run('/net/store/nbp/projects/phasesim/src_Arushi/matlab/addScriptPaths.m') 
      compl_fs_mask = load_untouch_nii('/net/store/nbp/projects/phasesim/workdir/Arushi/20150423gridjob/compl_fs_mask.nii');
      load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150423gridjob/new_tract_space.mat');
      first_voxel = ((this.params.ProbtrackX.split-1) * 1000) +1;
      probtrackx2Call = '. /usr/share/fsl/5.0/etc/fslconf/fsl.sh; /usr/share/fsl/5.0/bin/probtrackx2 -x /work/agarg/wimat2optionsinglemask/seedmask.nii -V 1 -l --onewaycondition --omatrix2 --target2=/work/agarg/wimat2optionsinglemask/termmask.nii -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --avoid=/net/store/nbp/projects/phasesim/databases/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/ca01_1_FA_thr0.12.nii.gz --stop=/work/agarg/wimat2optionsinglemask/termmask.nii --forcedir --opd -s /net/store/nbp/projects/phasesim/databases/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/bedpost.bedpostX/merged -m /net/store/nbp/projects/phasesim/databases/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/bedpost.bedpostX/nodif_brain_mask --dir=/work/agarg/wimat2optionsinglemask/data';      
      connmat = zeros(size(new_tract_space,1),size(new_tract_space,1));
      mkdir('/work/agarg/wimat2optionsinglemask')
      
      if this.params.ProbtrackX.split ~= 53
        for i=first_voxel:first_voxel+this.params.ProbtrackX.numberPerSplit-1
                %Fetching coordinates of the ith voxel
                x = new_tract_space(i,1:3);
  
                %Check if the ith voxel is 1
                if compl_fs_mask.img(x(1),x(2),x(3)) == 1
    
                    % set the voxel to 0 for term mask and save mask for probtrackx2 run
                    compl_fs_mask.img(x(1),x(2),x(3)) = 0;     
                    delete('/work/agarg/wimat2optionsinglemask/termmask.nii')
                    save_untouch_nii(compl_fs_mask, '/work/agarg/wimat2optionsinglemask/termmask');
                    
                    %set the value in the matrix back to 1 for the next
                    %runs
                    compl_fs_mask.img(x(1),x(2),x(3)) = 1;
                    seed_voxel = compl_fs_mask;
                    seed_voxel.img(:) = 0;
                    seed_voxel.img(x(1),x(2),x(3)) = 1;
                    delete('/work/agarg/wimat2optionsinglemask/seedmask.nii')
                    save_untouch_nii(seed_voxel, '/work/agarg/wimat2optionsinglemask/seedmask');
                    
                    % call for probtrackx2
                    status = system(probtrackx2Call);
                    
                    if status
                        %save error to file including mentioning the voxel
                        %number
                    else
                        fdt_matrix = load('/work/agarg/wimat2optionsinglemask/data/fdt_matrix2.dot');
                        
                        for j = 1 : size(fdt_matrix,1)
                            connmat(i,fdt_matrix(j,2)) = connmat(i,fdt_matrix(j,2)) + fdt_matrix(j,3);
                        end
                        rmdir('/work/agarg/wimat2optionsinglemask/data','s');
                    end
                             
                    
                   
                else
                    disp('voxel', i, 'is not populated in the fs mask');
                end
        end
      else
        for i = first_voxel:size(new_tract_space)
            
                %Fetching coordinates of the ith voxel
                x = new_tract_space(i,1:3);
  
                %Check if the ith voxel is 1
                if compl_fs_mask.img(x(1),x(2),x(3)) == 1
    
                   % set the voxel to 0 for term mask and save mask for probtrackx2 run
                    compl_fs_mask.img(x(1),x(2),x(3)) = 0;     
                    delete('/work/agarg/wimat2optionsinglemask/termmask.nii')
                    save_untouch_nii(compl_fs_mask, '/work/agarg/wimat2optionsinglemask/termmask');
                    
                    %set the value in the matrix back to 1 for the next
                    %runs
                    compl_fs_mask.img(x(1),x(2),x(3)) = 1;
                    
                    seed_voxel = compl_fs_mask;
                    seed_voxel.img(:) = 0;
                    seed_voxel.img(x(1),x(2),x(3)) = 1;
                    delete('/work/agarg/wimat2optionsinglemask/seedmask.nii')
                    save_untouch_nii(seed_voxel, '/work/agarg/wimat2optionsinglemask/seedmask');
                    
                   % call for probtrackx2
                    status = system(probtrackx2Call);
                    if status
                        %save error to file including mentioning the voxel
                        %number
                    else
                          fdt_matrix = load('/work/agarg/wimat2optionsinglemask/data/fdt_matrix2.dot');
                          waytotal = load('work/agarg/wimat2optionsinglemask/data/waytotal');
                          waytotaltotal = waytotaltotal + waytotal;
                          
                        
                        for j = 1 : size(fdt_matrix,1)
                            connmat(i,fdt_matrix(j,2)) = connmat(i,fdt_matrix(j,2)) + fdt_matrix(j,3);
                        end
                        rmdir('/work/agarg/wimat2optionsinglemask/data','s');
                    end
                    
                else
                    disp('voxel', i, 'is not populated in the fs mask');
                end
        end
      end
      
      save(['/net/store/projects/phasesim/workdir/Arushi/20150423gridjob/connmat' this.params.ProbtrackX.split], 'connmat', '-v7.3'); 
      save(['/net/store/projects/phasesim/workdir/Arushi/20150423gridjob/waytotal' this.params.ProbtrackX.split], 'waytotaltotal');
      rmdir('/work/agarg','s');
                    
      save(fullfile(this.workpath,[this.currJobid '.mat']));
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      for i = 1:53
          temp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150423gridjob/connmat' i '.mat']);
          connmat = connmat + temp;
          save('/net/store/nbp/projects/phasesim/workdir/Arushi/20150423gridjob/connmat' , 'connmat','-v7.3')
          
          tempwaytotal = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150423gridjob/waytotal' i '.mat']);
          waytotal = waytotal + tempwaytotal;
          save('/net/store/nbp/projects/phasesim/workdir/Arushi/20150423gridjob/waytotal','waytotal');
      end      
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

