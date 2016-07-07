function h = plotVoxelsInBrain(roi_coords_tract_space, subj_fs_mask_loc)

%% load nifti header to calculate transform matrix (from Tractography Space to Freesurfer Space):
if isempty(subj_fs_mask_loc)
    header = cbiReadNiftiHeader('/net/store/nbp/projects/phasesim/databases/SC_Bastian/mosaic_2015_02_18/compl_fs_mask.nii');
else 
    header = cbiReadNiftiHeader(subj_fs_mask_loc);
end
srow_mat = [header.srow_x; header.srow_y; header.srow_z];

%% load voxel coordinates in tractography space:
roi_coords_fs_space = [roi_coords_tract_space ones(size(roi_coords_tract_space,1),1)] * srow_mat';

%% plot:
h = plot3(roi_coords_fs_space(:,1),roi_coords_fs_space(:,2),roi_coords_fs_space(:,3),'o','MarkerFaceColor','b','markersize',3);

end