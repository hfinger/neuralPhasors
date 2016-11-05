classdef ProbtrackXallsubjects < Gridjob
  %ProbtrackX Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    
  end
  
  methods
    
    %% Subclass Constructor: initialize standard parameters:
    function this = ProbtrackXallsubjects(varargin)
      
      % Call superclass constructor:
      this = this@Gridjob(varargin{:});
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: define standard parameters for the job %%%%
      
      this.params.ProbtrackXallsubjects.split = num2cell(1:100);
      this.params.ProbtrackXallsubjects.splitPerSubject = 100;
      this.params.ProbtrackXallsubjects.subjectId = num2cell([1:4 6:13 15 17:22]);
      
      %%%% END EDIT HERE:                                          %%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      this = this.init(varargin{:});
      
    end
    
    %% Start: is executed before all individual parameter jobs are started
    function startJob(this)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%% START EDIT HERE: do some preprocessing (not parallel) %%%%
      
      disp('ProbtrackXall:this excecutes before the parallel job is started');
      
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
      
      datapaths = dataPaths();
      datapaths.workdir
      
      %call program to wait
      waitForServer();
      
      %system call & calculations for retrieving memory of work directory
      [status, result] = system('df -h 2>/dev/null | grep /work');
      disp(result);
      
      if status
        error(['df -h | grep /work;' 'result' result]);
      end
      
      resultwospace = strread(result, '%s');
      disp(resultwospace);
      freememory = resultwospace(end-2);
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
%       voxelcount = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20160210Allsubjecttracking/voxelcount_FA_thr_012.mat');
      
      %set flag to move to local if memory is more than 500 MB
      if strfind(powerfreememory, 'M')
        %check if it is more than 500MB
        if numfreememory > 500
          storelocalflag = 1;  %set flag to move to local
        else
          storelocalflag = 0;
        end
      elseif strfind(powerfreememory, 'G')
        storelocalflag = 1; %set flag to move to local
      elseif strfind(powerfreememory, 'T')
        storelocalflag = 1;
      else
        storelocalflag = 0;
      end
      
      % change for all subject tracking 20150727
      subjectId = this.params.ProbtrackXallsubjects.subjectId;
      if subjectId < 10;
        caNum = ['0' num2str(subjectId)];
      else
        caNum = num2str(subjectId);
      end
      
      %%%
      %       compl_fs_mask = load_untouch_nii([datapaths.workdir
      %       '/Arushi/20150423gridjob/compl_fs_mask.nii']);
      compl_fs_mask = load_untouch_nii([datapaths.workdir '/Arushi/20160210Allsubjecttracking/ca' caNum 'FA_masks_FA_thr_012compl_fs_mask.nii']);
      voxelIndexInImg = find(compl_fs_mask.img);
      %change end for all subject tracking 20150727
      
      voxelCount = size(voxelIndexInImg,1);      
      splitPerSubject = this.params.ProbtrackXallsubjects.splitPerSubject;
      
      currentSplit = this.params.ProbtrackXallsubjects.split;
      numberVoxPerSplit = ceil(voxelCount / splitPerSubject);
      if currentSplit == splitPerSubject
        % if we we are in the last split of this subject, there are less
        % voxels remaining to process:
%         numberVoxThisSplit = voxelCount - numberVoxPerSplit * (1 - splitPerSubject);
        numberVoxThisSplit = voxelCount - ((numberVoxPerSplit-1) * splitPerSubject) ;
      else
        numberVoxThisSplit = numberVoxPerSplit;
      end
      
      fdt_matrix = cell(numberVoxThisSplit,1);
      fdt_matrix_lengths = cell(numberVoxThisSplit,1);
      waytotal = zeros(numberVoxThisSplit,1);
      
      %change start for all subject tracking
      FAmask = load_untouch_nii([datapaths.databases '/Bastian_DTI/dti_data/ca' caNum '_1/ca' caNum '_1_FA.nii.gz']);
      FAmask.img(FAmask.img <= 0.1)=0;
      FAmask.img(FAmask.img~=0)=1;
      FAmask.img = ~FAmask.img;
      
      if storelocalflag
        [~, username] = system('whoami');
        userfolder = fullfile('/work', username(1:end-1));
        temppath = [userfolder '/probtrackx_tmp/' this.params.Gridjob.jobname '_jobid_' num2str(this.currJobid)];
        mkdir(temppath);
        save_untouch_nii(FAmask, [temppath '/ca' caNum 'FA_thr_0.1.nii']);
        copyfile([datapaths.databases '/Bastian_DTI/dti_data/ca' caNum '_1/bedpost.bedpostX']...
          ,[temppath '/bedpost.bedpostX']);
        capath = [temppath '/ca' caNum 'FA_thr_0.1.nii'];
        bedpostpath = [temppath '/bedpost.bedpostX'];
      else
        temppath = [this.workpath '/probtrackx_tmp/' this.params.Gridjob.jobname '_jobid_' num2str(this.currJobid)];
        if ~exist(temppath);
          mkdir(temppath);
        end
        save_untouch_nii(FAmask, [temppath '/ca' caNum 'FA_thr_0.1.nii']);
        capath = [temppath '/ca' caNum 'FA_thr_0.1.nii'];
        bedpostpath = [datapaths.databases '/Bastian_DTI/dti_data/ca' caNum '_1/bedpost.bedpostX'];
      end
      
      first_voxel = ((currentSplit - 1) * numberVoxPerSplit) + 1;
      for k=1:numberVoxThisSplit
        
        voxId = k+first_voxel-1;
        
        %Check if the ith voxel is 1
        if compl_fs_mask.img(voxelIndexInImg(voxId)) ~= 1
          error(['voxel', num2str(voxId), 'is not populated in the fs mask']);
        end
        
        % set the voxel to 0 for term mask and save mask for probtrackx2 run
        compl_fs_mask.img(voxelIndexInImg(voxId)) = 0;
        delete([temppath '/termmask.nii'])
        save_untouch_nii(compl_fs_mask, [temppath '/termmask']);
        
        %set the value in the matrix back to 1 for the next
        %runs
        compl_fs_mask.img(voxelIndexInImg(voxId)) = 1;
        seed_voxel = compl_fs_mask;
        seed_voxel.img(:) = 0;
        seed_voxel.img(voxelIndexInImg(voxId)) = 1; %last change 20150728
        delete([temppath '/seedmask.nii'])
        save_untouch_nii(seed_voxel, [temppath '/seedmask']);
        
        probtrackx2Call = ['. /usr/share/fsl/5.0/etc/fslconf/fsl.sh; /usr/share/fsl/5.0/bin/probtrackx2'...
          ' -x ' temppath '/seedmask.nii'...
          ' -V 1'...
          ' -l --onewaycondition --omatrix2'...
          ' --target2=' temppath '/termmask.nii'...
          ' -c 0.2 -S 2000 --steplength=0.5'...
          ' -P 5000'...
          ' --rseed=12345'...
          ' --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 '...
          ' --avoid=' capath...
          ' --stop=' temppath '/termmask.nii --forcedir --opd'...
          ' -s ' bedpostpath '/merged'...
          ' --ompl'...
          ' -m ' bedpostpath '/nodif_brain_mask'...
          ' --dir=' temppath '/data' num2str(k) ];
        % call for probtrackx2
        system('free');
        disp(datetime);
        status = system(probtrackx2Call);
        disp(['status:' num2str(status)]);
        system('free');
        disp(datetime);
        
        if status
          disp(['error in subjnum' num2str(subjectId) 'voxel:' num2str(voxId)]);
          ME = MException('MyComponent:ProbtrackXrun', ...
            'Error in subjnum %s voxel %s run %s',num2str(subjectId), num2str(voxId), num2str(k));
          throw(ME)
        else
          fdt_matrix{k} = load([temppath '/data' num2str(k) '/fdt_matrix2.dot']);
          fdt_matrix_lengths{k} = load([temppath '/data' num2str(k) '/fdt_matrix2_lengths.dot']);
          
          % correct the target indices due to removing the seed voxel from the target mask:
          indicesToIncrease = find(fdt_matrix{k}(:,2) >= k);
          fdt_matrix{k}(indicesToIncrease,2) = fdt_matrix{k}(indicesToIncrease,2) + 1;
          fdt_matrix_lengths{k}(indicesToIncrease,2) = fdt_matrix_lengths{k}(indicesToIncrease,2) + 1;
          
          fdt_matrix{k}(:,1) = k;
          fdt_matrix_lengths{k}(:,1) = k;
          
          waytotal(k) = load([temppath '/data' num2str(k) '/waytotal']);
          
          %                             for j = 1 : size(fdt_matrix,1)
          %                                 connmat(k,fdt_matrix(j,2)) = connmat(k,fdt_matrix(j,2)) + fdt_matrix(j,3);
          %                             end
          rmdir([temppath '/data' num2str(k)],'s');
        end
        
      end
      
      fdt_matrix = cat(1,fdt_matrix{:});
      connmat = spconvert(fdt_matrix); %#ok<NASGU>
      
      fdt_matrix_lengths = cat(1,fdt_matrix_lengths{:});
      distmat = spconvert(fdt_matrix_lengths); %#ok<NASGU>
      
      save([this.workpath '/connmat_subj_' caNum '_split_' num2str(currentSplit)], 'connmat', '-v7.3');
      save([this.workpath '/distmat_subj_' caNum '_split_' num2str(currentSplit)], 'distmat', '-v7.3');
      save([this.workpath '/waytotal_subj_' caNum '_split_' num2str(currentSplit)], 'waytotal');
      rmdir(temppath,'s');
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
      %         minIdx = (i-1)*this.params.ProbtrackX.numberVoxThisSplit+1;
      %         maxIdx = i*this.params.ProbtrackX.numberVoxThisSplit;
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

