function graclusNcut(recursiveSplit, clusterCount, subjRange, InputPath, OutputPath, normBy)
if recursiveSplit
    
    for subjNum = subjRange
        
        if normBy == 'sum'
                    compSimMat = load([InputPath '/compSimMat/' 'compSimMatsum' num2str(subjNum)]);
        else

        compSimMat = load([InputPath '/compSimMat/' 'compSimMat' num2str(subjNum)]);
        end
        compSimMat = compSimMat.expEucMat;
        
        InterimOutputPath = [OutputPath '/graclusIntermedResult/gracIntResultnormBy' normBy num2str(subjNum) ];
        if ~exist(InterimOutputPath, 'dir')
            mkdir(InterimOutputPath)
        end
        
        [clusterIdPerVoxel, clusterIdPerVoxelCurrent, largestClusterId, cutValue] = applyClustering( compSimMat, clusterCount, InterimOutputPath );
        
        stepWiseClusters = clusterIdPerVoxel;
        allClusters = clusterIdPerVoxelCurrent;
        
        FinalOutputPath = [OutputPath 'normby' normBy];
        
        
        if ~exist(FinalOutputPath, 'dir')
            mkdir(FinalOutputPath);
        end
        
        
        save([FinalOutputPath '/graclusResultnormBy' normBy num2str(subjNum)], 'stepWiseClusters', 'allClusters', 'largestClusterId', 'cutValue');
        
    end
else
    search_steps = 80;
    [clusterIdPerVoxel, cutValue] = graclus(high_res_tracts, clusterCount, 0, search_steps, 0);
    largestClusterId = [];
end
end