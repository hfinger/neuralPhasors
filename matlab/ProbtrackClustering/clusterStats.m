function clusterStats()

for subjNum = 18:22
    if subjNum < 10
        caNum = ['0' num2str(subjNum)];
    else
        caNum = num2str(subjNum);
    end
    maskPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/ca' caNum 'FA_masks_FA_thr_012compl_fs_mask.nii'];
    useVoxelPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/' num2str(subjNum) '/'];
    if exist(maskPath, 'file')
    mask = load_untouch_nii(maskPath);
    else
        continue;
    end
    Allvoxels = find(mask.img);
    
    useVoxelNum = 1;
    while 1
        
        useVoxelFinalPath = [useVoxelPath 'useVoxelIdx' num2str(useVoxelNum) '.mat'];
        if ~exist(useVoxelFinalPath, 'file')
            break;
        end
        useVoxelIdx = load(useVoxelFinalPath);
        useVoxelIdx = useVoxelIdx.useVoxelIdx;
        Allvoxels = Allvoxels(useVoxelIdx);
        
        useVoxelNum = useVoxelNum + 1;
    end
        
        
    for mode = 1:4
        switch mode
            
            case 1
                modeText = 'fullconn';
                startId = 2;
                
            case 2
                modeText = 'fullcos';
                startId = 2;
                
            case 3
                modeText = 'fsconn';
                startId = 67;
                
            case 4
                modeText = 'fscos';
                startId = 67;
                
        end
        if subjNum == 18
            if ~strcmp(modeText,'fscos')
                continue;
            end
        end
        clusterCoMAll = cell(1,1000);
        clusterSizeAll = cell(1,1000);
        clusterSecondMomentAll = cell(1,1000);
        splitClusterId = cell(1,1000);
        splitClusterIdBySize = cell(1,1000);
        savepath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150920StatsAndMetrics/' num2str(subjNum) '/' modeText];
        
        for clusterNum = startId:1001
            loadfilepath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/' num2str(subjNum) '/' modeText '/' modeText 'Cluster' num2str(clusterNum)];
            disp(['subj' num2str(subjNum) modeText num2str(clusterNum)]);
            %    clusterConnmat = load([filepath '/' 'clusterConnmat' num2str(i) '.mat']);
            %    clusterConnmat = clusterConnmat.clusterConnmat;
            %    clusterConnmat = bsxfun(@rdivide, clusterConnmat, sum(clusterConnmat,2));
            %    clusterConnmat = 1./(clusterConnmat);
            %    betCentr = betweenness_wei(clusterConnmat);
            
            Cluster = load([loadfilepath '.mat']);
            Cluster = Cluster.Cluster;
            
            clusterVoxels = arrayfun(@(x) find(Cluster==x), 1:clusterNum, 'UniformOutput', false);
            clusterCoords = cellfun(@(x) Allvoxels(x,:), clusterVoxels, 'UniformOutput', false);
            clusterSize = cellfun(@(x) length(x), clusterVoxels, 'UniformOutput', false);
            clusterCoMAll{clusterNum} = cellfun(@(x) mean(x,1), clusterCoords, 'UniformOutput', false);
            
            if ~isempty(clusterCoMAll{clusterNum-1})
                lastCoords = cell2mat(clusterCoMAll{clusterNum-1}(1:end)');
                currentCoords = cell2mat(clusterCoMAll{clusterNum}(1:end-1)');
                
                splitClusterId{clusterNum} = find(sum(lastCoords~=currentCoords,2));
            end
            
            if ~isempty(clusterSizeAll{clusterNum-1})
                [~, splitClusterIdBySize{clusterNum}] = max(clusterSizeAll{clusterNum-1});
            end
            
            clusterSecondMoment = cellfun(@(x) mean(sqrt(sum(bsxfun(@minus,x,mean(x,1)).^2,2))), clusterCoords, 'UniformOutput', false);
            
            clusterSizeAll{clusterNum} = cell2mat(clusterSize);
            clusterSecondMomentAll{clusterNum} = cell2mat(clusterSecondMoment);
            
        end
        
        mkdir(savepath)
        save(fullfile(savepath, 'stats.mat'), 'clusterSizeAll', 'clusterSecondMomentAll','clusterCoMAll','splitClusterId','splitClusterIdBySize');
    end
end
end