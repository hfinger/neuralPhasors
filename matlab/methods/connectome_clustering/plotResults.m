close all;
clear all;
clc;

path_prefix = '/net/store/nbp/projects/phasesim/databases/DTI_subject01/';
% load_path = '/net/store/nbp/projects/phasesim/workdir/Holger/20160330_GraclusClustering/new_results_weight0.5_decay0.25/clusters.mat';
% voxel_path = [path_prefix 'tractographyData/voxel_coords.mat'];


load_path = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/decay-1weigh0.5connnormbysumthresh100/graclusResultnormBysumthresh100subj1.mat';
voxel_path = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/decay-1weigh0.5/FinalCoord/FinalCoord1.mat';

%% load clustering data:
clusters = load(load_path);
stepWiseClusters = clusters.stepWiseClusters;
largestClusterId = clusters.largestClusterId;
cutValue = clusters.cutValue;

%% load surface data:

freesurfer_roi_ids = load([path_prefix 'tractographyData/freesurfer_roi_ids.mat']);
freesurfer_roi_ids = freesurfer_roi_ids.freesurfer_roi_ids;

voxel_coords = load(voxel_path);
voxel_coords = voxel_coords.FinalCoord;
voxel_count = size(voxel_coords, 1);

header = cbiReadNiftiHeader([path_prefix 'surfaceFsData/compl_fs_mask.nii']);
g = gifti([path_prefix 'surfaceFsData/ca01_1_structcortex_8196.surf.gii']);

%% load our surface
% 
% % g_ours_lh = gifti('/net/store/nbp/projects/phasesim/databases/Bastian_DTI/new/freesurfer/test_ca01_1/surf/lh.white.native.gii');
% g_ours_lh = gifti('/net/store/nbp/projects/phasesim/databases/Bastian_DTI/new/ca01_1/lh.white.gii');
% g_ours_rh = gifti('/net/store/nbp/projects/phasesim/databases/Bastian_DTI/new/ca01_1/rh.white.gii');
% 
% % apply transformations:
% surfaceTrafo = g_ours_lh.mat;
% vertices = g_ours_lh.vertices;
% vertices = (surfaceTrafo*[vertices ones(size(vertices,1),1)]')';
% g_ours_lh.vertices = vertices(:,1:3);
% g_ours_lh.mat = eye(4,4);
% 
% % construct full brain surface:
% % g_ours = gifti();
% % g_ours.faces = g_ours_lh.faces;
% % g_ours.mat = eye(4,4);
% % g_ours.vertices = g_ours_lh.vertices;
% 
% % g_ours_lh.vertices(:,1) = g_ours_lh.vertices(:,1) + 13;
% % g_ours_lh.vertices(:,2) = g_ours_lh.vertices(:,2) - 20;
% 
% % plot full brain:
% az=270;
% el=30;
% figure(1);
% clf;
% hg_ours = plot(g_ours_lh);
% alpha(hg_ours,0.1);
% set(gca,'xtick',[],'ytick',[])
% hold on;
% hh = plotVoxelsInBrain(header, voxel_coords);
% % h = plot3(voxel_coords(:,1),voxel_coords(:,2),voxel_coords(:,3),'.','markersize',3);
% hold off;
% view(az,el);

%% loop over iterations:

az=0;
el=25;
figure(2);
clf;

for iteration = 1:1000
  data = stepWiseClusters(:,iteration);
  clf;
  
  title(['iter ' num2str(iteration) ' cutVal ' num2str(cutValue(iteration))])
  
  % load brain surface and plot it:
  hg = plot(g);
  alpha(hg,0.05);
  set(gca,'xtick',[],'ytick',[])
  
  % plot some ROI's:
  hold on;
  hh = plotVoxelsInBrain(header, voxel_coords(data==iteration,:));
  hh = plotVoxelsInBrain(header, voxel_coords(data==largestClusterId(iteration),:));
  hold off
  
  for t=1:150
    az=az+2;
    view(az,el);
    drawnow;
        pause(0.01);
  end
  
end


%% plot single iteration
iteration = 1;
figure(1);
data = stepWiseClusters(:,iteration);
clf;
hg = plot(g);
alpha(hg,0.1);
set(gca,'xtick',[],'ytick',[])
% plot some ROI's:
hold on;
for i=1:randi(1000,50)
  hh = plotVoxelsInBrain(header, voxel_coords(data==i,:));
end
hold off


%% analyze distribution of roi sizes:

iteration=500;
tmp=stepWiseClusters(:,iteration);
a = unique(tmp);
out = [a,histc(tmp(:),a)];
hist(double(out(:,2)),[0:200])
disp(['minimum cluster size: ' num2str(min(out(:,2)))])

%%
figure(3)
plot(cutValue)

%%
largestClusterId(1)=1; %fix first value
numVoxInLargestCluster = zeros(999,1);
numVoxInNewCluster1 = zeros(999,1);
numVoxInNewCluster2 = zeros(999,1);
for iteration=1:999
  clusterId = largestClusterId(iteration);
  numVoxInLargestCluster(iteration) = sum(stepWiseClusters(:,iteration) == clusterId);
  numVoxInNewCluster1(iteration) = sum(stepWiseClusters(:,iteration+1) == clusterId);
  numVoxInNewCluster2(iteration) = sum(stepWiseClusters(:,iteration+1) == iteration+1);
end
sizeRatio = abs(numVoxInNewCluster1-numVoxInNewCluster2)./(numVoxInNewCluster1+numVoxInNewCluster2);

%%
figure(4)
clf;
plot(sizeRatio)
hold on;
plot(cutValue)

%%
scatter(sizeRatio(801:900),cutValue(800:899))