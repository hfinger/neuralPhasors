close all;
clear all;
clc;

savepath_prefix = '/net/store/nbp/projects/phasesim/workdir/Holger/20160330_GraclusClustering/';

%% load tractography data:
disp('load tractography data')
path_prefix = '/net/store/nbp/projects/phasesim/databases/DTI_subject01/';

freesurfer_roi_ids = load([path_prefix 'tractographyData/freesurfer_roi_ids.mat']);
freesurfer_roi_ids = freesurfer_roi_ids.freesurfer_roi_ids;

voxel_coords = load([path_prefix 'tractographyData/voxel_coords.mat']);
voxel_coords = voxel_coords.voxel_coords;
voxel_count = size(voxel_coords, 1);

header = cbiReadNiftiHeader([path_prefix 'surfaceFsData/compl_fs_mask.nii']);
g = gifti([path_prefix 'surfaceFsData/ca01_1_structcortex_8196.surf.gii']);

high_res_tracts = load([path_prefix 'tractographyData/high_res_tracts.mat']);
high_res_tracts = high_res_tracts.high_res_tracts;

distance_matrix = load('/net/store/nbp/projects/phasesim/src_enes2/code&data/distances/250_distance_matrix');
distance_matrix = distance_matrix.distance_matrix;

%% make the matrix symmetrical
disp('make the matrix symmetrical')
high_res_tracts = bsxfun(@rdivide,high_res_tracts,sum(high_res_tracts,2));
high_res_tracts = high_res_tracts + high_res_tracts';

%% convert to cosine similarity matrix:
useCosineSimilarity = true;
if useCosineSimilarity
  disp('convert to cosine similarity matrix')
  tmp = high_res_tracts*high_res_tracts';
  tmp = bsxfun(@rdivide,tmp, sqrt(sum(high_res_tracts.^2,2)) );
  high_res_tracts = bsxfun(@rdivide,tmp, sqrt(sum(high_res_tracts.^2,1)) );
  clear tmp;
end

%% make main diagonal zero
disp('make main diagonal zero')
high_res_tracts(logical(eye(size(high_res_tracts)))) = 0;

%% calculate exponential decay for the distance matrix
disp('calculate exponential decay for the distance matrix')
decay_constant = 0.25;
distance_matrix = spfun(@exp, -distance_matrix * decay_constant);

%% normalize the matrices:
disp('normalize the matrices')
distance_matrix = bsxfun(@rdivide,distance_matrix,sum(distance_matrix,2));
high_res_tracts = bsxfun(@rdivide, high_res_tracts, sum(high_res_tracts,2));

%% combine DTI data with distance matrix:
disp('combine DTI data with distance matrix')
distance_weight = 0;
connectivity_weight = 1 - distance_weight;
high_res_tracts = (distance_weight*distance_matrix) + (connectivity_weight*high_res_tracts);

%% do clustering
disp('do clustering')
recursive_splitting = true;
cluster_count = 1000;
if recursive_splitting
  [clusterIdPerVoxel, largestClusterId, cutValue] = applyClustering( high_res_tracts, cluster_count );
else
  search_steps = 80;
  [clusterIdPerVoxel, cutValue] = graclus(high_res_tracts, cluster_count, 0, search_steps, 0);
  largestClusterId = [];
end

%% save
disp('save')

results_dir = ['weight' num2str(distance_weight) '_decay' num2str(decay_constant)];
if useCosineSimilarity
  results_dir = ['cosineSim_' results_dir];
end
if ~recursive_splitting
  results_dir = ['notRecursive_' results_dir];
end
results_dir = ['results_' results_dir];
results_path = fullfile(savepath_prefix, results_dir);
mkdir(results_path)
save(fullfile(results_path,'clusters.mat'), 'clusterIdPerVoxel', 'largestClusterId', 'cutValue')
save(fullfile(results_path,'high_res_tracts.mat'), 'high_res_tracts', '-v7.3')

%% plot results
if 1
  iteration = 1000;
  figure(1);
  data = clusterIdPerVoxel(:,iteration);
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
end

%% plot results
if 1
  az=0;
  el=25;
  figure(2);
  
  for iteration = 1
    data = clusterIdPerVoxel(:,iteration);
    clf;
    
    title(num2str(iteration))
    
    % load brain surface and plot it:
    hg = plot(g);
    alpha(hg,0.1);
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
      %     pause(0.05);
    end
    
  end
end
