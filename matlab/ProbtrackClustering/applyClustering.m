function [clusterIdPerVoxel, clusterIdPerVoxelCurrent, largestClusterId, cutValue, threshFactor] = applyClustering( high_res_tracts, target_cluster_count, recursiveSplit, threshFactor )
% log
% 20160410 - agarg - copied from matlab/methods/connectome_clustering
%                  - edited return arguments to include clusterIdPerVoxelCurrent
%                  - edited return arguments to include OutputPath to save
%                    intermediate results


%% use graclus for clustering
% G: Graph adjacency matrix
% k: number of clusters
% cutType: 0 for NCut, 1 for RAssoc (default is NCut)
% l: number olf local search steps (default is 0)
% spectral: 1 for spectral clustering at coarsest level

split_count = 2;
cut_type = 0;
search_steps = 80;
spectral = 0;

clusterIdPerVoxel = zeros(size(high_res_tracts,1),target_cluster_count,'uint16');
largestClusterId = zeros(target_cluster_count,1);
cutValue = zeros(target_cluster_count,1);

% init first iter (all voxels belong to cluster with id=1):
clusterIdPerVoxelCurrent = ones(size(high_res_tracts,1),1,'uint16');
clusterIdPerVoxel(:,1) = clusterIdPerVoxelCurrent;

%% iterative clustering
if recursiveSplit
    for i = 2:target_cluster_count
        disp(['number clusters: ' num2str(i)]);
        
        % find largest cluster:
        largestClusterId(i) = mode(clusterIdPerVoxelCurrent);
        voxelIdsInLargestCluster = find(clusterIdPerVoxelCurrent == largestClusterId(i));
         clusterSC = high_res_tracts(voxelIdsInLargestCluster, voxelIdsInLargestCluster);
        
       
                    [partition, cutValue(i)] = graclus(clusterSC, split_count, cut_type, search_steps, spectral);
               
                    % give one of the two partitions a new cluster id (which is the current iteration):
        clusterIdPerVoxelCurrent(voxelIdsInLargestCluster(partition==1)) = i;
        
        % save into full matrix:
        clusterIdPerVoxel(:,i) = clusterIdPerVoxelCurrent;
             
        
    end
    
    %% Non-recursive splitting
    
else
    search_steps = 80;
    [clusterIdPerVoxel, cutValue] = graclus(high_res_tracts, target_cluster_count, 0, search_steps, 0);
    largestClusterId = [];
end

end

