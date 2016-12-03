function [number_of_outliers, number_of_connectedComponents] = connectedComp(voxel_coords)
% This function calculates the number of connected components and number of
% outlier voxels

% [numOutliersVoxels, NumObjects] = evaluate(partition, voxel_coords, clusterID)

allowed_gap_size = 3;

cluster_voxel_coords = voxel_coords;
clusterCenter = mean(cluster_voxel_coords,1);
%clusterSecondMoment(cluster_ID) = mean(sqrt(sum(bsxfun(@minus,cluster_voxel_coords,clusterCenter).^2,2)).^2);

img3d = zeros(max(voxel_coords,[],1));
cluster_voxel_indices = sub2ind(size(img3d),cluster_voxel_coords(:,1),cluster_voxel_coords(:,2),cluster_voxel_coords(:,3));
img3d(cluster_voxel_indices) = 1;
if allowed_gap_size>0
  img3d = imdilate(img3d,ones(allowed_gap_size*2+1,allowed_gap_size*2+1,allowed_gap_size*2+1));
end
cc = bwconncomp(img3d,26);

number_of_outliers = 0;
for c=1:length(cc.PixelIdxList)
  component_indicies = intersect(cc.PixelIdxList{c},cluster_voxel_indices);
  if c>1
    number_of_outliers = number_of_outliers + length(component_indicies);
  end
end

number_of_connectedComponents = cc.NumObjects;
disp(['number connected components: ' num2str(number_of_connectedComponents)])
disp(['number outlier voxels: ' num2str(number_of_outliers)])
end

