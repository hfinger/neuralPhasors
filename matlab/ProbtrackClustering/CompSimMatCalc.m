function CompSimMatCalc(subjRange, InputPath, OutputPath, decayParam, WeighingFactor, calcDistMat, normBy, onlyChangeNorm)
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
%% Calculate normalised and reduced connmat
for subjNum = subjRange
    if ~onlyChangeNorm
    % Load connectivity matrix
    S = load([InputPath '/connMat/' 'connmatSubj' num2str(subjNum) '.mat']);
    S = S.connmat;
    
    % Find indices where rows sum to zero, columns sum to zero and both
    zeroRowVoxelIdx    = find(sum(S,2)==0);
    zeroColumnVoxelIdx = find(sum(S,1)==0);
    zeroBothVoxelIdx   = intersect(zeroRowVoxelIdx, zeroColumnVoxelIdx);
    
    % Remove rows and columns only where both sum to zero
    S(zeroBothVoxelIdx,:) = [];
    S(:,zeroBothVoxelIdx) = [];
    
    if ~exist([OutputPath '/zeroVoxelIdx/'], 'dir')
        mkdir([OutputPath '/zeroVoxelIdx/']);
    end
    save([OutputPath '/zeroVoxelIdx/' 'zeroRowVoxelIdx' num2str(subjNum)],'zeroRowVoxelIdx');
    save([OutputPath '/zeroVoxelIdx/' 'zeroColumnVoxelIdx' num2str(subjNum)],'zeroColumnVoxelIdx');
    save([OutputPath '/zeroVoxelIdx/' 'zeroBothVoxelIdx' num2str(subjNum)], 'zeroBothVoxelIdx');
    clear zeroRowVoxelIdx;
    clear zeroColumnVoxelIdx;
    clear zeroBothVoxelIdx;
    
    % Normalise by sum of rows
    S = bsxfun(@rdivide,S,sum(S,2));
    
    % Convert all NaNs to zero
    S(isnan(S)) = 0;
    
    % Make symmetric
    S = S + S';
    % Make diagonal zero
    S = S-diag(diag(S));
    
    if ~exist([OutputPath '/SSymDiag0/'], 'dir')
        mkdir([OutputPath '/SSymDiag0/']);
    end
    save([OutputPath '/SSymDiag0/' 'SSymDiag0subj' num2str(subjNum)], 'S', '-v7.3');
    
    else
    
        S = load([OutputPath '/SSymDiag0/' 'SSymDiag0subj' num2str(subjNum)]);
        S = S.S;
    end
    
    switch normBy
        case 'mean'
            % Normalise by mean of whole matrix;
            S = S/mean(S(:));
            
        case 'sum'
            S = S/sum(S(:)) * size(S,1);
    end
    
    if ~exist([OutputPath '/SNormandRed/'], 'dir')
        mkdir([OutputPath '/SNormandRed/']);
    end
    save([OutputPath '/SNormandRed/' 'SNormandRed' normBy num2str(subjNum)], 'S', '-v7.3');
    clear S;
    
end
%% Calculate Euclidean distance matrix for old subject1 data

for subjNum = subjRange
    
    if ~onlyChangeNorm
    if subjNum < 10;
        caNum = ['0' num2str(subjNum)];
    else
        caNum = num2str(subjNum);
    end
    
    % Load complete grey matter masks
    complFSMaskfilename = ['ca' caNum 'FA_masks_FA_thr_012compl_fs_mask'];
    complFSMask = load_untouch_nii([InputPath '/complFSMask/' complFSMaskfilename '.nii']);
    
    % Load indices where both rows and columns sum to zero
    zeroBothVoxelIdx = load([OutputPath '/zeroVoxelIdx/' 'zeroBothVoxelIdx' num2str(subjNum)]);
    zeroBothVoxelIdx = zeroBothVoxelIdx.zeroBothVoxelIdx;
    
    % find all grey matter voxel indices and coordinates
    compInd = find(complFSMask.img);
    [compSubx, compSuby, compSubz] = ind2sub(size(complFSMask.img), compInd);
    FinalCoord = [compSubx compSuby compSubz];
    FinalCoord(zeroBothVoxelIdx,:) = [];
    
    % Save
    if ~exist([OutputPath '/FinalCoord/'], 'dir')
        mkdir([OutputPath '/FinalCoord/']);
    end
    save([OutputPath '/FinalCoord/' 'FinalCoord' num2str(subjNum)], 'FinalCoord');
    
    totVoxel = 1:size(FinalCoord,1);
    distMat = sparse(zeros(size(FinalCoord,1)));
    clear compSubx compSuby compSubz compInd compSub;
    clear complFSMask;
    
    if calcDistMat
        % Calculate Euclidean distance matrix
        for n = totVoxel
            tmp = sqrt(bsxfun(@minus, FinalCoord(totVoxel,1),FinalCoord(n,1)).^2 + ...
                bsxfun(@minus, FinalCoord(totVoxel,2),FinalCoord(n,2)).^2 + ...
                bsxfun(@minus, FinalCoord(totVoxel,3),FinalCoord(n,3)).^2 );
            
            keep = find(tmp<=6);
            
            distMat(n, keep) = tmp(keep);
            
        end
        if ~exist([OutputPath '/distMat/'], 'dir')
            mkdir([OutputPath '/distMat/']);
        end
        save([OutputPath '/distMat/' 'distMat' num2str(subjNum)], 'distMat', '-v7.3');
    else
        load([OutputPath '/distMat/' 'distMat' num2str(subjNum)]);
    end
    expEucMat = spfun(@exp, -distMat * decayParam);
    clear distMat;
    
    % Make diagonal zero
    expEucMat = expEucMat-diag(diag(expEucMat));
    
    save([OutputPath '/expEucMat/' 'expEucMat' num2str(subjNum)], 'expEucMat', '-v7.3');
    else
        expEucMat = load([OutputPath '/expEucMat/' 'expEucMat' num2str(subjNum)]);
        expEucMat = expEucMat.expEucMat;
    end
    switch normBy
        case 'mean'
            % Normalise by mean of whole matrix;
            expEucMat = expEucMat/mean(expEucMat(:));
            
        case 'sum'
            expEucMat = expEucMat/sum(expEucMat(:)) * size(expEucMat,1);
    end
    if ~exist([OutputPath '/expEucMat/'], 'dir')
        mkdir([OutputPath '/expEucMat/']);
    end
    save([OutputPath '/expEucMat/' 'normexpEucMat' normBy num2str(subjNum)], 'expEucMat', '-v7.3');
    
    % Add Spatial Distance and Connectivity matrix
    load([OutputPath '/SNormandRed/' 'SNormandRed' normBy num2str(subjNum) '.mat']);
    
    expEucMat = expEucMat*WeighingFactor + S*(1-WeighingFactor);
    
    
    if ~exist([OutputPath '/compSimMat/'], 'dir')
        mkdir([OutputPath '/compSimMat/']);
    end
    save([OutputPath '/compSimMat/' 'compSimMat' normBy num2str(subjNum)],'expEucMat','-v7.3');
    clear expEucMat;
    clear S;
    clear compInd
    clear compSub
    clear n
    clear totVoxel
end
end