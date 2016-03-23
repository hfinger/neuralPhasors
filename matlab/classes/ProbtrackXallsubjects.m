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
            %change 20150727
            this.params.ProbtrackXallsubjects.split = num2cell(116); %change
            this.params.ProbtrackXallsubjects.numSubjects = 22;
            this.params.ProbtrackXallsubjects.splitPerSubject = 50;
            %       this.params.ProbtrackX.numberPerSplit = 1000;  %change
            %       this.params.ProbtrackX.numberRegions = 52442;  %change
            %change end 20150727
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
            
            [~, username] = system('whoami');
            userfolder = fullfile('/work', username(1:end-1));
            jobdirectory = fullfile(userfolder, 'tempProbtrackXallsubjects');
            temppath = fullfile(jobdirectory, num2str(this.currJobid));
            mkdir(temppath);
            
            datapaths = dataPaths();
            datapaths.workdir
            
            
            
            
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
            voxelcount = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/voxelcount_FA_thr_012.mat');
            
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
            
            % change for all subject tracking 20150727
            subjectNum =  this.params.ProbtrackXallsubjects.split/this.params.ProbtrackXallsubjects.splitPerSubject;
            if floor(subjectNum) ~= subjectNum
                subjectNum = ceil(subjectNum);
            else
                
            end
            
            if subjectNum < 10;
                caNum = ['0' num2str(subjectNum)];
            else
                caNum = num2str(subjectNum);
            end
            
            if ~(subjectNum == 5 | subjectNum == 14 | subjectNum == 16)
            %%%
            %       compl_fs_mask = load_untouch_nii([datapaths.workdir
            %       '/Arushi/20150423gridjob/compl_fs_mask.nii']);
            compl_fs_mask = load_untouch_nii([datapaths.workdir '/Arushi/20150727Allsubjecttracking/ca' caNum 'FA_masks_FA_thr_012compl_fs_mask.nii']);
            voxelIndex = find(compl_fs_mask.img);
            %change end for all subject tracking 20150727
            
            %change start for all subject tracking 20150728
            FAmask = load_untouch_nii([datapaths.databases '/Bastian_DTI/dti_data/ca' caNum '_1/ca' caNum '_1_FA.nii.gz']);
            FAmask.img(FAmask.img <= 0.1)=0;
            FAmask.img(FAmask.img~=0)=1;
            FAmask.img = ~FAmask.img;
            save_untouch_nii(FAmask, [temppath '/ca' caNum 'FA_thr_0.1.nii']);
            
            voxelCount = voxelcount.voxelcount_FA_thr_012(subjectNum);
            numberPerSplit = ceil(voxelCount/50);
            connmat = zeros(numberPerSplit,voxelCount);
            waytotal = zeros(numberPerSplit,1);
            
            
            
            %       load([datapaths.workdir '/Arushi/20150423gridjob/new_tract_space.mat']);
            %       connmat = zeros(this.params.ProbtrackX.numberPerSplit,size(new_tract_space,1));
            %       waytotal = zeros(this.params.ProbtrackX.numberPerSplit,1);
            
            if storelocalflag
                %             copyfile([datapaths.databases '/SC_Bastian/mosaic_2015_02_18/ca01_1_dti/ca01_1_FA_thr0.12.nii.gz']...
                
                copyfile([datapaths.databases '/Bastian_DTI/dti_data/ca' caNum '_1/bedpost.bedpostX']...
                    ,[temppath '/bedpost.bedpostX']);
                capath = [temppath '/ca' caNum 'FA_thr_0.1.nii'];
                
                bedpostpath = [temppath '/bedpost.bedpostX'];
            else
                capath = [datapaths.databases '/Bastian_DTI/dti_data/ca' caNum '_1/ca' caNum '_1_FA.nii.gz'];
                bedpostpath = [datapaths.databases '/Bastian_DTI/dti_data/ca' caNum '_1/bedpost.bedpostX'];
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
            subjPartition = this.params.ProbtrackXallsubjects.split - ((subjectNum-1)*this.params.ProbtrackXallsubjects.splitPerSubject);
            first_voxel = ((subjPartition - 1)*numberPerSplit) +1;
            for k=1:numberPerSplit
                
                
                i = k+first_voxel-1;
                
                % break if i is greater than
                if i > voxelCount
                    break;
                end
                
               
                
                %Check if the ith voxel is 1
                if compl_fs_mask.img(voxelIndex(i)) == 1
                    
                    % set the voxel to 0 for term mask and save mask for probtrackx2 run
                    compl_fs_mask.img(voxelIndex(i)) = 0;
                    delete([temppath '/termmask.nii'])
                    save_untouch_nii(compl_fs_mask, [temppath '/termmask']);
                    
                    %set the value in the matrix back to 1 for the next
                    %runs
                    compl_fs_mask.img(voxelIndex(i)) = 1;
                    seed_voxel = compl_fs_mask;
                    seed_voxel.img(:) = 0;
                    seed_voxel.img(voxelIndex(i)) = 1; %last change 20150728
                    delete([temppath '/seedmask.nii'])
                    save_untouch_nii(seed_voxel, [temppath '/seedmask']);
                    
                    % call for probtrackx2
                    status = system(probtrackx2Call);
                    disp(['status:' num2str(status)]);
                    
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
            connmat = sparse(connmat);
            
            save([this.workpath '/connmat' num2str(this.params.ProbtrackXallsubjects.split)], 'connmat', '-v7.3');
            save([this.workpath '/waytotal' num2str(this.params.ProbtrackXallsubjects.split)], 'waytotal');
            rmdir(temppath,'s');
            rmdir(jobdirectory);
            status = rmdir(userfolder);
            
            %%%% END EDIT HERE:                                %%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
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

