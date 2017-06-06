% To calculate Composite Similarity Matrices
%20160407 - created
%20160408 - -commenting;
%           -added to make diagonal zero in the final matrix;
%           -remove line to make symmetric before removing sum zeros and use
%           it later after normalisation by sum of rows;
%           -remove line to remove voxels which row sum to zero
%           -add line to save list of indices which row sum to zero or
%           column sum to zero
% 20160412  - Make diagonal zero in respective matrices before normalising
%           by mean
%           - shift exponential conversion of euclidean matrix outside the
%           for loop and save distance matrix before euclidean conversion
%           to use later; include parameter to say when to calc distMat
%           - include parameter to choose normalisation
% 20160415 - Change order of processing. Throw exception if diagonal is
%            nonzero
%          - Include Cosine Similarity parameter and processing
% 20160504 - Changed to a script to make it easier instead of having too
%            many input parameters
%          - added inputdlg

%% Parameters
prompt = {'Enter Weighing Factor range:',...
    'Enter subject range:',...
    'Enter Decay Parameter:',...
    'Do you want to use Cosine Similarity("true" or "false")?',...
    'Do you want to calculate S matrix("true" or "false")?',...
    'Do you want to calculate distance matrix("true" or "false")?',...
    'Enter individual normalisation factor("sum" or "mean")',...
    'Do you only want to change normalisation("true" or "false")?',...
    'Enter whole normalisation("WholeMax")',...
    'Enter Path to pick up data:',...
    'Enter Main Path to save data:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'(0:0.1:0.9)',...
    '[1:4, 6:13, 15, 17:22]',...
    '-1', ...
    'false',...
    'true',...
    'true',...
    'sum',...
    'true',...
    'WholeMax',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/20160417_Allconnmat',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

%% Processing
if ~isempty(answer)
    subjRange = str2num(answer{2});
    decayParam = str2num(answer{3});
    useCosineSim = str2num(answer{4});
    calcSMat = str2num(answer{5});
    calcDistMat = str2num(answer{6});
    normBy = answer{7};
    onlyChangeNorm = str2num(answer{8});
    WholeNormText = answer{9};
    InputPath = answer{10};
    MainPath = answer{11};
    WeighFacRange = str2num(answer{1});
    if useCosineSim
        cosText = 'cos';
    else
        cosText = 'conn';
    end
    CosPath = [MainPath '/' cosText];
    %% Calculate normalised and reduced connmat
    if calcSMat
        for subjNum = subjRange
            if ~onlyChangeNorm
                % Load connectivity matrix
                S = load([InputPath '/connmatSubj' num2str(subjNum) '.mat']);
                S = S.connmat;
                
                % Check if diagonal zero
                diagonal = diag(S);
                DiagNonZero = find(diagonal);
                clear diagonal;
                
                if DiagNonZero
                    disp(['error in subjnum' num2str(subjNum)]);
                    ME = MException('MyComponent:ProbtrackXrun', ...
                        'Error in subjnum %s',num2str(subjectId));
                    throw(ME)
                end
                clear DiagNonZero;
                
                % Find indices where rows sum to zero, columns sum to zero and both
                zeroRowVoxelIdx    = find(sum(S,2)==0);
                zeroColumnVoxelIdx = find(sum(S,1)==0);
                zeroBothVoxelIdx   = intersect(zeroRowVoxelIdx, zeroColumnVoxelIdx);
                clear zeroRowVoxelIdx;
                clear zeroColumnVoxelIdx;
                
                % Remove rows and columns only where both sum to zero
                S(zeroBothVoxelIdx,:) = [];
                S(:,zeroBothVoxelIdx) = [];
                clear zeroBothVoxelIdx;
                
                if ~exist([MainPath '/zeroVoxelIdx/'], 'dir')
                    mkdir([MainPath '/zeroVoxelIdx/']);
                end
                save([MainPath '/zeroVoxelIdx/' 'zeroRowVoxelIdx' num2str(subjNum)],'zeroRowVoxelIdx');
                save([MainPath '/zeroVoxelIdx/' 'zeroColumnVoxelIdx' num2str(subjNum)],'zeroColumnVoxelIdx');
                save([MainPath '/zeroVoxelIdx/' 'zeroBothVoxelIdx' num2str(subjNum)], 'zeroBothVoxelIdx');
                clear zeroRowVoxelIdx;
                clear zeroColumnVoxelIdx;
                clear zeroBothVoxelIdx;
                
                % Normalise by sum of rows
                S = bsxfun(@rdivide,S,sum(S,2));
                
                % Convert all NaNs to zero
                S(isnan(S)) = 0;
                
                % Make symmetric
                S = S + S';
                
                
                if ~exist([MainPath '/SSymDiag0/'], 'dir')
                    mkdir([MainPath '/SSymDiag0/']);
                end
                save([MainPath '/SSymDiag0/' 'SSymDiag0subj' num2str(subjNum)], 'S', '-v7.3');
                
                if useCosineSim
                    tmp = S*S';
                    tmp = bsxfun(@rdivide,tmp, sqrt(sum(S.^2,2)) );
                    S = bsxfun(@rdivide,tmp, sqrt(sum(S.^2,1)) );
                    clear tmp;
                    
                    if ~exist([CosPath '/ScosSim/'], 'dir')
                        mkdir([CosPath '/ScosSim/']);(0)
                    end
                    save([CosPath '/ScosSim/' 'ScosSimSubj' num2str(subjNum)], 'S', '-v7.3');
                end
            else
                if useCosineSim
                    S = load([MainPath '/SSymDiag0/' 'SSymDiag0subj' num2str(subjNum)]);
                    S = S.S;
                    tmp = S*S';
                    tmp = bsxfun(@rdivide,tmp, sqrt(sum(S.^2,2)) );
                    S = bsxfun(@rdivide,tmp, sqrt(sum(S.^2,1)) );
                    clear tmp;
                    
                    if ~exist([CosPath '/ScosSim/'], 'dir')
                        mkdir([CosPath '/ScosSim/']);
                    end
                    save([CosPath '/ScosSim/' 'ScosSimSubj' num2str(subjNum)], 'S', '-v7.3');
                    
                else
                    S = load([MainPath '/SSymDiag0/' 'SSymDiag0subj' num2str(subjNum)]);
                    S = S.S;
                end
            end
            
            switch normBy
                case 'mean'
                    % Normalise by mean of whole matrix;
                    S = S/mean(S(:));
                    
                case 'sum'
                    S = S/sum(S(:)) * size(S,1);
            end
            
            if ~exist([CosPath '/SNormandRed/'], 'dir')
                mkdir([CosPath '/SNormandRed/']);
            end
            save([CosPath '/SNormandRed/' 'SNormandRed' normBy num2str(subjNum)], 'S', '-v7.3');
            clear S;
            
        end
    end
    
    %% Calculate Euclidean distance matrix
    if calcDistMat
        
        for subjNum = subjRange
            
            if ~onlyChangeNorm
                if subjNum < 10;
                    caNum = ['0' num2str(subjNum)];
                else
                    caNum = num2str(subjNum);
                end
                
                % Load complete grey matter masks
                complFSMaskfilename = ['ca' caNum 'FA_masks_FA_thr_012compl_fs_mask'];
                complFSMask = load_untouch_nii(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSub/complFSMask/' complFSMaskfilename '.nii']);
                clear caNum;
                
                % Load indices where both rows and columns sum to zero
                zeroBothVoxelIdx = load([MainPath '/zeroVoxelIdx/' 'zeroBothVoxelIdx' num2str(subjNum)]);
                zeroBothVoxelIdx = zeroBothVoxelIdx.zeroBothVoxelIdx;
                
                % find all grey matter voxel indices and coordinates
                compInd = find(complFSMask.img);
                [compSubx, compSuby, compSubz] = ind2sub(size(complFSMask.img), compInd);
                FinalCoord = [compSubx compSuby compSubz];
                FinalCoord(zeroBothVoxelIdx,:) = [];
                clear zeroBothVoxelIdx;
                clear compInd;
                clear compSubx;
                clear compSuby;
                clear compSubz;
                
                
                % Save
                if ~exist([MainPath '/FinalCoord/'], 'dir')
                    mkdir([MainPath '/FinalCoord/']);
                end
                save([MainPath '/FinalCoord/' 'FinalCoord' num2str(subjNum)], 'FinalCoord');
                
                totVoxel = 1:size(FinalCoord,1);
                distMat = sparse(zeros(size(FinalCoord,1)));
                clear compSubx compSuby compSubz compInd compSub;
                clear complFSMask;
                
                % Calculate Euclidean distance matrix
                for n = totVoxel
                    tmp = sqrt(bsxfun(@minus, FinalCoord(totVoxel,1),FinalCoord(n,1)).^2 + ...
                        bsxfun(@minus, FinalCoord(totVoxel,2),FinalCoord(n,2)).^2 + ...
                        bsxfun(@minus, FinalCoord(totVoxel,3),FinalCoord(n,3)).^2 );
                    
                    keep = find(tmp<=6);
                    
                    distMat(n, keep) = tmp(keep);
                    
                end
                clear FinalCoord;
                clear totVoxel;
                clear tmp;
                clear n;
                if ~exist([MainPath '/distMat/'], 'dir')
                    mkdir([MainPath '/distMat/']);
                end
                save([MainPath '/distMat/' 'distMat' num2str(subjNum)], 'distMat', '-v7.3');
                expEucMat = spfun(@exp, -distMat * decayParam);
                clear distMat;
                
                % Make diagonal zero
                expEucMat = expEucMat-diag(diag(expEucMat));
                
                
                if ~exist([MainPath '/expEucMat/'], 'dir')
                    mkdir([MainPath '/expEucMat/']);
                end
                
                save([MainPath '/expEucMat/' 'expEucMatdecay' num2str(decayParam) 'subj' num2str(subjNum)], 'expEucMat', '-v7.3');
            else
                expEucMat = load([MainPath '/expEucMat/' 'expEucMatdecay' num2str(decayParam) 'subj' num2str(subjNum)]);
                expEucMat = expEucMat.expEucMat;
            end
            switch normBy
                case 'mean'
                    % Normalise by mean of whole matrix;
                    expEucMat = expEucMat/mean(expEucMat(:));
                    
                case 'sum'
                    expEucMat = expEucMat/sum(expEucMat(:)) * size(expEucMat,1);
            end
            if ~exist([MainPath '/expEucMat/'], 'dir')
                mkdir([MainPath '/expEucMat/']);
            end
            save([MainPath '/expEucMat/' 'normexpEucMat' normBy 'decay' num2str(decayParam) 'subj' num2str(subjNum)], 'expEucMat', '-v7.3');
            clear expEucMat;
            
            
        end
    end
    %% Calculate composite similarity matrix
    for subjNum = subjRange
        
        load([MainPath '/expEucMat/' 'normexpEucMat' normBy 'decay' num2str(decayParam) 'subj' num2str(subjNum)]);
        
        load([CosPath '/SNormandRed/' 'SNormandRed' normBy num2str(subjNum)]);
        
        for WeighingFactor = WeighFacRange
            
            OutputPath = [MainPath '/' cosText '/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor)]  ;
            
            
            % Add Spatial Distance and Connectivity matrix
            
            compSimMat = expEucMat*WeighingFactor + S*(1-WeighingFactor);
            
            
            if ~exist([OutputPath '/compSimMat/'], 'dir')
                mkdir([OutputPath '/compSimMat/']);
            end
            save([OutputPath '/compSimMat/' 'compSimMat' normBy num2str(subjNum)],'compSimMat','-v7.3');
            
            if strcmp(WholeNormText, 'WholeMax')
                compSimMat = compSimMat/max(compSimMat(:));
            end
            
            if ~exist([OutputPath '/compSimMat' WholeNormText '/'], 'dir')
                mkdir([OutputPath '/compSimMat' WholeNormText '/']);
            end
            save([OutputPath '/compSimMat' WholeNormText '/' 'compSimMat' WholeNormText  normBy num2str(subjNum)],'compSimMat','-v7.3');
            clear compSimMat;
        end
        
        clear expEucMat;
        clear S;
    end
end