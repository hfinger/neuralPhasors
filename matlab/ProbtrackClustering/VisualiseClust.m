%VisualiseClust

%% load tractography data:
cluster_roi_ids = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fullcos/fullcosCluster2.mat');
cluster_roi_ids = cluster_roi_ids.Cluster;

voxel_coords = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20160209_High_Res_SC/tractographyData/voxel_coords.mat');
voxel_coords = voxel_coords.voxel_coords;

subj_surf_file_loc = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/ca01_1_structcortex_8196.surf.gii';
subj_fs_mask_loc = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/mosaic_2015_02_18/compl_fs_mask.nii';

g = testPlot(subj_surf_file_loc, subj_fs_mask_loc, cluster_roi_ids, voxel_coords);