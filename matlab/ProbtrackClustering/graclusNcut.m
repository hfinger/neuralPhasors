function graclusNcut(recursiveSplit, clusterCount, subjRange, InputPath, OutputPath, normBy, threshFactor)
    
    for subjNum = subjRange
        
        compSimMat = load([InputPath '/compSimMatWholeMax/' 'compSimMatWholeMax' normBy num2str(subjNum)]);
       
        compSimMat = compSimMat.compSimMat;
        
        compSimMat = compSimMat*(threshFactor);
        
              
        [clusterIdPerVoxel, clusterIdPerVoxelCurrent, largestClusterId, cutValue]...
            = applyClustering( compSimMat, clusterCount, recursiveSplit );
        
        stepWiseClusters = clusterIdPerVoxel;
        allClusters = clusterIdPerVoxelCurrent;
        
       
        FinalOutputPath = [OutputPath 'normby' normBy 'thresh' num2str(threshFactor) ];
        
        
        if ~exist(FinalOutputPath, 'dir')
            mkdir(FinalOutputPath);
        end
        
        
        save([FinalOutputPath '/graclusResultnormBy' normBy 'thresh' num2str(threshFactor) 'clust' num2str(clusterCount) ...
            'subj' num2str(subjNum)], 'stepWiseClusters', 'allClusters', 'largestClusterId', 'cutValue');
        
    end

end