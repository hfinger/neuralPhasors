function ret = calcCommonClusters(voxelIndByCluster1, voxelIndByCluster2)

IntersectVal = cell(length(find(~cellfun(@isempty, voxelIndByCluster1))), length(find(~cellfun(@isempty, voxelIndByCluster2))));
IntersectInd1 = cell(size(IntersectVal));
IntersectInd2 = cell(size(IntersectVal));
IntersectCount = zeros(size(IntersectVal));
IntersectRatio1 = zeros(size(IntersectVal));
IntersectRatio2 = zeros(size(IntersectVal));
ret.MaxIntRatio1 = zeros(length(find(~cellfun(@isempty, voxelIndByCluster1))),1);
ret.MaxIntRatio2 = zeros(length(find(~cellfun(@isempty, voxelIndByCluster2))),1);
% IndMaxIntRatio1 = cell(size(ret.MaxIntRatio1));
% IndMaxIntRatio2 = cell(size(ret.MaxIntRatio2));

for clustNum1 = 1:length(find(~cellfun(@isempty, voxelIndByCluster1)))
    for clustNum2 = 1:length(find(~cellfun(@isempty, voxelIndByCluster2)))
        
        [IntersectVal{clustNum1,clustNum2}, IntersectInd1{clustNum1, clustNum2}, IntersectInd2{clustNum1,clustNum2}]...
            = intersect(voxelIndByCluster1{clustNum1}, voxelIndByCluster2{clustNum2});
        
        IntersectCount(clustNum1, clustNum2) = length(IntersectInd1{clustNum1,clustNum2});
        
        IntersectRatio1(clustNum1, clustNum2) = IntersectCount(clustNum1, clustNum2)/length(voxelIndByCluster1{clustNum1});
        IntersectRatio2(clustNum1, clustNum2) = IntersectCount(clustNum1, clustNum2)/length(voxelIndByCluster2{clustNum2});
  
            
        
    end
end
clear IntersectVal;
clear IntersectCount;
clear IntersectInd1;
clear IntersectInd2
for clustNum1 = 1:length(find(~cellfun(@isempty, voxelIndByCluster1)));
   ret.MaxIntRatio1(clustNum1) = max(IntersectRatio1(clustNum1,:));
   
%     IndMaxIntRatio1{clustNum1} = find(IntersectRatio1(clustNum1,:) == ret.MaxIntRatio1(clustNum1));
end
for clustNum2 =  1:length(find(~cellfun(@isempty, voxelIndByCluster2)))
        ret.MaxIntRatio2(clustNum2) = max(IntersectRatio2(:,clustNum2));
        
%         IndMaxIntRatio2{clustNum2} = find(ret.IntersectRatio2(:,clustNum2) == ret.MaxIntRatio2(clustNum2));
end

end
    

