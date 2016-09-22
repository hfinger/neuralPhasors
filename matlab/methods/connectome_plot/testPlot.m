%% load tractography data:

path_prefix = '/net/store/nbp/projects/phasesim/databases/DTI_subject01/';

freesurfer_roi_ids = load([path_prefix 'tractographyData/freesurfer_roi_ids.mat']);
freesurfer_roi_ids = freesurfer_roi_ids.freesurfer_roi_ids;

voxel_coords = load([path_prefix 'tractographyData/voxel_coords.mat']);
voxel_coords = voxel_coords.voxel_coords;

header = cbiReadNiftiHeader([path_prefix 'surfaceFsData/compl_fs_mask.nii']);
g = gifti([path_prefix 'surfaceFsData/ca01_1_structcortex_8196.surf.gii']);


%% create new figure and plot the surface and some ROI's:
figure(1); 
clf;

% load brain surface and plot it:
hg = plot(g);
alpha(hg,0.1);
set(gca,'xtick',[],'ytick',[])

% plot some ROI's:
hold on;
h1 = plotVoxelsInBrain(header, voxel_coords(freesurfer_roi_ids==1,:));
h2 = plotVoxelsInBrain(header, voxel_coords(freesurfer_roi_ids==2,:));
h3 = plotVoxelsInBrain(header, voxel_coords(freesurfer_roi_ids==3,:));
hold off