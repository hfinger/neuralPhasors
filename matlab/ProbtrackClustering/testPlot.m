function g = testPlot(subj_surf_file_loc, subj_fs_mask_loc, cluster_roi_ids, voxel_coords)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% include libraries:
addpath('/net/store/nbp/projects/phasesim/workdir/Arushi/20160209_High_Res_SC/include/NIfTImatlab/matlab')
addpath('/net/store/nbp/projects/phasesim/workdir/Arushi/20160209_High_Res_SC/include/gifti-1.4')



%% create new figure and plot the surface and some ROI's:
figure(1); 
clf;

% load brain surface and plot it:
if isempty(subj_surf_file_loc)
g = gifti('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/ca01_1_structcortex_8196.surf.gii');
else
    g = gifti(subj_surf_file_loc);
end
hg = plot(g);
alpha(hg,0.1);
set(gca,'xtick',[],'ytick',[])

% plot some ROI's:
hold on;
h1 = plotVoxelsInBrain(voxel_coords(cluster_roi_ids==1,:), subj_fs_mask_loc);
h2 = plotVoxelsInBrain(voxel_coords(cluster_roi_ids==2,:), subj_fs_mask_loc);
hold off
end

