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
      
      [~, username] = system('whoami');
      userfolder = fullfile('/work', username(1:end-1));
      jobdirectory = fullfile(userfolder, 'tempProbtrackX');
      temppath = fullfile(jobdirectory, num2str(this.currJobid));
      mkdir(temppath);
      
      datapaths = dataPaths();
      datapaths.workdir
      
      
      first_voxel = ((this.params.ProbtrackX.split-1) * 1000) +1;
     
      
      %call program to wait
      waitForServer();
      
      %system call & calculations for retrieving memory of work directory
      [status, result] = system('df -h | grep /work'); 
            disp(result);
            
            resultwospace = strread(result, '%s');
            disp(resultwospace);
            freememory = resultwospace(4);
            freememory = cell2mat(freememory);
            

%       resultwospace = textscan(result, '%s');
%       disp(resultwospace);
%       
%       freememory = resultwospace{1,1};
%       freememory = freememory(4);
%       freememory = cell2mat(freememory);
      
      disp(freememory);
      numfreememory = freememory(1:regexp(freememory,'\D')-1);
      powerfreememory = freememory(regexp(freememory,'\D'));
      disp(powerfreememory);
      disp(numfreememory);
      
           
        %set flag to move to local if memory is more than 500 MB
      if powerfreememory == 'M'
          %check if it is more than 500MB
          if numfreememory > 500
          storelocalflag = 1;  %set flag to move to local
          end
      elseif powerfreememory == 'G'
          storelocalflag = 1; %set flag to move to local
      elseif powerfreememory == 'T'
          storelocalflag = 1;
      else
          storelocalflag = 0;
      end
       
      compl_fs_mask = load_untouch_nii([datapaths.workdir '/Arushi/20150423gridjob/compl_fs_mask.nii']);
      load([datapaths.workdir '/Arushi/20150423gridjob/new_tract_space.mat']);
      connmat = zeros(this.params.ProbtrackX.numberPerSplit,size(new_tract_space,1));
      waytotal = zeros(this.params.ProbtrackX.numberPerSplit,1);
    
          if storelocalflag 
            copyfile([datapaths.databases '/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/ca01_1_FA_thr0.12.nii.gz']...
          , [temppath '/ca01_1_FA_thr0.12.nii.gz']);
      copyfile([datapaths.databases '/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/bedpost.bedpostX']...
          ,[temppath '/bedpost.bedpostX']);
            capath = [temppath '/ca01_1_FA_thr0.12.nii.gz'];
            bedpostpath = [temppath '/bedpost.bedpostX'];
          else
              capath = [datapaths.databases '/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/ca01_1_FA_thr0.12.nii.gz'];
              bedpostpath = [datapaths.databases '/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/bedpost.bedpostX'];
          end
          
           probtrackx2Call = ['. /usr/share/fsl/5.0/etc/fslconf/fsl.sh; /usr/share/fsl/5.0/bin/probtrackx2'...
        ' -x ' temppath '/seedmask.nii'...
        ' -V 1'...
        ' -l --onewaycondition --omatrix2'...
        ' --target2=' temppath '/termmask.nii'...
        ' -c 0.2 -S 2000 --steplength=0.5'...
        ' -P 5000'...
        ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 '...
        ' --avoid=' capath...
        ' --stop=' temppath '/termmask.nii --forcedir --opd'...
        ' -s ' bedpostpath '/merged'...
        ' -m ' bedpostpath '/nodif_brain_mask'...
        ' --dir=' temppath '/data'];
      
      for k=1:this.params.ProbtrackX.numberPerSplit
        
        i = k+first_voxel-1;
        
        % break if i is greater than
        if i > this.params.ProbtrackX.numberRegions
          break;
        end
        
        %Fetching coordinates of the ith voxel
        x = new_tract_space(i,1:3);
        
        %Check if the ith voxel is 1
        if compl_fs_mask.img(x(1),x(2),x(3)) == 1
          
          % set the voxel to 0 for term mask and save mask for probtrackx2 run
          compl_fs_mask.img(x(1),x(2),x(3)) = 0;
          delete([temppath '/termmask.nii'])
          save_untouch_nii(compl_fs_mask, [temppath '/termmask']);
          
          %set the value in the matrix back to 1 for the next
          %runs
          compl_fs_mask.img(x(1),x(2),x(3)) = 1;
          seed_voxel = compl_fs_mask;
          seed_voxel.img(:) = 0;
          seed_voxel.img(x(1),x(2),x(3)) = 1;
          delete([temppath '/seedmask.nii'])
          save_untouch_nii(seed_voxel, [temppath '/seedmask']);
          
          % call for probtrackx2
          status = system(probtrackx2Call);
          
          if status
            %save error to file including mentioning the voxel
            %number
          else
            fdt_matrix = load([temppath '/data/fdt_matrix2.dot']);
            waytotal(k) = load([temppath '/data/waytotal']);
            
            for j = 1 : size(fdt_matrix,1)
              connmat(k,fdt_matrix(j,2)) = connmat(k,fdt_matrix(j,2)) + fdt_matrix(j,3);
            end
            rmdir([temppath '/data'],'s');
          end
          
        else
          disp('voxel', i, 'is not populated in the fs mask');
        end
      end
      
      save([this.workpath '/connmat' num2str(this.params.ProbtrackX.split)], 'connmat', '-v7.3');
      save([this.workpath '/waytotal' num2str(this.params.ProbtrackX.split)], 'waytotal');
      rmdir(temppath,'s');
      rmdir(jobdirectory);
      status = rmdir(userfolder);
      
      %%%% END EDIT HERE:                                %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    end
    
    %% finishJobs: is executed once (after all individual parameter jobs are finished)
    function finishJobs(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some clean up and saving %%%%
      
      % TODO: do this in another job, because it requires much more RAM...
      
%       connmat = zeros(this.params.ProbtrackX.numberRegions,this.params.ProbtrackX.numberRegions);
%       waytotal = zeros(this.params.ProbtrackX.numberRegions,1);
%       
%       for i = 1:53
%         minIdx = (i-1)*this.params.ProbtrackX.numberPerSplit+1;
%         maxIdx = i*this.params.ProbtrackX.numberPerSplit;
%         
%         temp = load([this.workpath '/connmat' i '.mat']);
%         connmat(minIdx:maxIdx,:) = connmat(minIdx:maxIdx,:) + temp.connmat;
%         save([this.workpath '/connmat'] , 'connmat','-v7.3')
%         
%         tempwaytotal = load([this.workpath '/waytotal' i '.mat']);
%         waytotal(minIdx:maxIdx) = waytotal(minIdx:maxIdx) + tempwaytotal.waytotal;
%         save([this.workpath '/waytotal'],'waytotal');
%       end
      
      %%%% END EDIT HERE:                               %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Execute clean up of superclass:
      finishJobs@Gridjob(this)
      
    end
    
  end
  
end

