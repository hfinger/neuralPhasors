function calcCommonClusters(voxelIndByCluster1, voxelIndByCluster2, OutputPath, clusterNum)

IntersectVal = cell(length(find(~cellfun(@isempty, voxelIndByCluster1))), length(find(~cellfun(@isempty, voxelIndByCluster2))));
IntersectInd1 = cell(size(IntersectVal));
IntersectInd2 = cell(size(IntersectVal));
IntersectCount = zeros(size(IntersectVal));
IntersectRatio1 = zeros(size(IntersectVal));
IntersectRatio2 = zeros(size(IntersectVal));
MaxIntRatio1 = zeros(length(find(~cellfun(@isempty, voxelIndByCluster1))),1);
MaxIntRatio2 = zeros(length(find(~cellfun(@isempty, voxelIndByCluster2))),1);
IndMaxIntRatio1 = zeros(size(MaxIntRatio1));
IndMaxIntRatio2 = zeros(size(MaxIntRatio2));

for clustNum1 = 1:length(find(~cellfun(@isempty, voxelIndByCluster1)))
    for clustNum2 = 1:length(find(~cellfun(@isempty, voxelIndByCluster2)))
        
        [IntersectVal{clustNum1,clustNum2}, IntersectInd1{clustNum1, clustNum2}, IntersectInd2{clustNum1,clustNum2}]...
            = intersect(voxelIndByCluster1{clustNum1}, voxelIndByCluster2{clustNum2});
        
        IntersectCount(clustNum1, clustNum2) = length(IntersectInd1{clustNum1,clustNum2});
        
        IntersectRatio1(clustNum1, clustNum2) = IntersectCount(clustNum1, clustNum2)/length(voxelIndByCluster1{clustNum1});
        IntersectRatio2(clustNum1, clustNum2) = IntersectCount(clustNum1, clustNum2)/length(voxelIndByCluster2{clustNum2});
  
            
        
    end
end
for clustNum1 = 1:length(find(~cellfun(@isempty, voxelIndByCluster1)));
   MaxIntRatio1(clustNum1) = max(IntersectRatio1(clustNum1,:));
    IndMaxIntRatio1(clustNum1) = find(IntersectRatio1(clustNum1,:) == MaxIntRatio1(clustNum1));
end
for clustNum2 =  1:length(find(~cellfun(@isempty, voxelIndByCluster2)))
        MaxIntRatio2(clustNum2) = max(IntersectRatio2(:,clustNum2));
        
       
        IndMaxIntRatio2(clustNum2) = find(IntersectRatio2(:,clustNum2) == MaxIntRatio2(clustNum2));
end

if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end
save([OutputPath '/intersectionCluster' num2str(clusterNum)],...
    'IntersectVal', 'IntersectInd1', 'IntersectInd2', 'IntersectCount', 'IntersectRatio1', 'IntersectRatio2',...
    'MaxIntRatio1', 'MaxIntRatio2', 'IndMaxIntRatio1', 'IndMaxIntRatio2');
end
    

