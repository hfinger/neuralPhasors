clear all;

dpaths = dataPaths();
matlab_path = dpaths.sge_pathToAddScriptPaths;

%% load tractography data:
freesurfer_roi_ids = load(fullfile(matlab_path, 'methods/connectome/brain_plot_data/tractographyData/freesurfer_roi_ids.mat'));
freesurfer_roi_ids = freesurfer_roi_ids.freesurfer_roi_ids;

voxel_coords = load(fullfile(matlab_path, 'methods/connectome/brain_plot_data/tractographyData/voxel_coords.mat'));
voxel_coords = voxel_coords.voxel_coords;

%% calculate center of each roi:
roi_center_coords = zeros(66,3);
for k=1:66
    roi_center_coords(k,:) = mean(voxel_coords(freesurfer_roi_ids==k,:),1);
end

%% load nifti header to calculate transform matrix (from Tractography Space to Freesurfer Space):
header = cbiReadNiftiHeader(fullfile(matlab_path, 'methods/connectome/brain_plot_data/surfaceFsData/compl_fs_mask.nii'));
srow_mat = [header.srow_x; header.srow_y; header.srow_z];

%% load voxel coordinates in tractography space:
roi_coords_fs_space = [roi_center_coords ones(size(roi_center_coords,1),1)] * srow_mat';

%% create new figure and plot the surface and some ROI's:
figure(1); 
clf;

% load brain surface and plot it:
g = gifti(fullfile(matlab_path, 'methods/connectome/brain_plot_data/surfaceFsData/ca01_1_structcortex_8196.surf.gii'));
hg = plot(g);
alpha(hg,0.1);
set(gca,'xtick',[],'ytick',[])

%%
hold on;

h = plot3(roi_coords_fs_space(:,1),roi_coords_fs_space(:,2),roi_coords_fs_space(:,3),'o','MarkerFaceColor','b','markersize',3);

hold off;

%%

if true
    % plot some ROI's Voxels:
    hold on;
    h1 = plotVoxelsInBrain(voxel_coords(freesurfer_roi_ids==1,:));
    h2 = plotVoxelsInBrain(voxel_coords(freesurfer_roi_ids==2,:));
    h3 = plotVoxelsInBrain(voxel_coords(freesurfer_roi_ids==3,:));
    hold off
end