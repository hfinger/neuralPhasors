function h = plotVoxelsInBrain(roi_coords_tract_space)

%% load nifti header to calculate transform matrix (from Tractography Space to Freesurfer Space):
header = cbiReadNiftiHeader('surfaceFsData/compl_fs_mask.nii');
srow_mat = [header.srow_x; header.srow_y; header.srow_z];

%% load voxel coordinates in tractography space:
roi_coords_fs_space = [roi_coords_tract_space ones(size(roi_coords_tract_space,1),1)] * srow_mat';

%% plot:
h = plot3(roi_coords_fs_space(:,1),roi_coords_fs_space(:,2),roi_coords_fs_space(:,3),'o','MarkerFaceColor','b','markersize',3);

end