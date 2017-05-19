function GlobalCutValCalc( clustRange, compSimMatPath, OutputPath, clusterIdPath, thresholdFactor  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 
compSimMat = load(compSimMatPath);
compSimMat = compSimMat.compSimMat;

% renormalize so that 
renormFactor = sum(compSimMat(:))/size(compSimMat,1);
compSimMat = compSimMat / renormFactor;

%compSimMat = thresholdFactor * compSimMat;
%compSimMat = round(1000*compSimMat); % According to Graclus internal
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
    
    disp(num2str(clusterCount));
    
end
    
cumulGlobalCutVal = sum(clustGlobalCutVal);

if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end

save([OutputPath, 'GlobalCutVal'], 'cumulGlobalCutVal');
end