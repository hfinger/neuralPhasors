for clusterNum = 2:1000
    a.newClust{clusterNum} = stepWiseClusters(:,clusterNum) - stepWiseClusters(:,clusterNum-1);
    a.newClustInd{clusterNum} = find(a.newClust{clusterNum} ~=0);
    
    a.OldClustNum(clusterNum) = stepWiseClusters(a.newClustInd{clusterNum}(1),(clusterNum-1));
    a.NewClustNum(clusterNum) = stepWiseClusters(a.newClustInd{clusterNum}(1), clusterNum);
    
    a.newClustSizeOld(clusterNum) = length(voxelIndByCluster{a.OldClustNum(clusterNum), clusterNum});
    a.newClustSizeNew(clusterNum) = length(voxelIndByCluster{a.NewClustNum(clusterNum), clusterNum});
end

a.OldClustSize = a.newClustSizeOld + a.newClustSizeNew;
a.RatioSize = a.newClustSizeOld./a.newClustSizeNew;
for i = 2:1000 
if a.RatioSize(i) > 1 
a.RatioSize(i) = 1/a.RatioSize(i); 
end 
end