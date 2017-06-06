function h = plotVoxelsInBrain(roi_coords_tract_space, subj_fs_mask_loc)

%% load nifti header to calculate transform matrix (from Tractography Space to Freesurfer Space):
if isempty(subj_fs_mask_loc)
    header = cbiReadNiftiHeader('/net/store/nbp/projects/phasesim/databases/SC_Bastian/mosaic_2015_02_18/compl_fs_mask.nii');
else 
    header = load_untouch_header_only('/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSub/complFSMask/ca01FA_masks_FA_thr_012compl_fs_mask.nii');
end
srow_mat = [header.hist.srow_x; header.hist.srow_y; header.hist.srow_z];

%% load voxel coordinates in tractography space:
roi_coords_fs_space = [roi_coords_tract_space ones(size(roi_coords_tract_space,1),1)] * srow_mat';

%% plot:
h = plot3(roi_coords_fs_space(:,1),roi_coords_fs_space(:,2),roi_coords_fs_space(:,3),'o','MarkerFaceColor','b','markersize',3);

end