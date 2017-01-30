%
       voxelsInClust = (ClusterIdsforVoxel(:,clusterCount) == clustNum);
        tmp = compSimMat(voxelsInClust, ~voxelsInClust);
        links = sum(tmp(:));
        degree = sum(compSimMat(voxelsInClust,:),:);