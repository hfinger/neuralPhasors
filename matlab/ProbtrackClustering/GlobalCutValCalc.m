function GlobalCutValCalc( clustRange, compSimMatPath, GraclusPath, clusterIdPath, thresholdFactor  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
compSimMat = load(compSimMatPath);
compSimMat = compSimMat.compSimMat;
compSimMat = thresholdFactor * compSimMat;
compSimMat = round(1000*compSimMat); % According to Graclus internal
ClusterIdsforVoxel = load(clusterIdPath);
ClusterIdsforVoxel = ClusterIdsforVoxel.stepWiseClusters;
clustGlobalCutVal = zeros(size(ClusterIdsforVoxel,2));

for clusterCount = clustRange
    
    clustGlobalCutVal(:,clusterCount) = clustGlobalCutVal(:,clusterCount-1);
    
    newClust = ClusterIdsforVoxel(:,clusterCount) - ClusterIdsforVoxel(:,clusterCount-1);
    newClustInd = find(newClust ~=0);
    
    OldClustNum = ClusterIdsforVoxel(newClustInd(1),(clusterCount-1));
    NewClustNum = ClusterIdsforVoxel(newClustInd(1), clusterCount);
    
    voxelsInOldClust = (ClusterIdsforVoxel(:,clusterCount) == OldClustNum);
    tmp = compSimMat(voxelsInOldClust, ~voxelsInOldClust);
    links = sum(tmp(:));
    tmp = compSimMat(voxelsInOldClust, :);
    degree = sum(tmp(:));
    clustGlobalCutVal(OldClustNum,clusterCount) = links/degree;
    
    
    voxelsInNewClust = (ClusterIdsforVoxel(:,clusterCount) == NewClustNum);
     tmp = compSimMat(voxelsInNewClust, ~voxelsInNewClust);
    links = sum(tmp(:));
    tmp = compSimMat(voxelsInNewClust, :);
    degree = sum(tmp(:));
    clustGlobalCutVal(NewClustNum,clusterCount) = links/degree;
    
end
    
cumulGlobalCutVal = sum(clustGlobalCutVal);

save([GraclusPath '/globalCutVal'], 'clustGlobalCutVal', 'cumulGlobalCutVal');
end

%   for clustNum = 1:clusterCount
%         
%         voxelsInClust = (ClusterIdsforVoxel(:,clusterCount) == clustNum);
%         tmp = compSimMat(voxelsInClust, ~voxelsInClust);
%         links = sum(tmp(:));
%         tmp = compSimMat(voxelsInClust, :);
%         degree = sum(tmp(:));
%         
%         globalCutVal(clusterCount) = globalCutVal(clusterCount) + links/degree;
%         
%     end
