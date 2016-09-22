close all;
clear all;
clc;

results_dir = '/net/store/nbp/projects/phasesim/workdir/Holger/20160330_GraclusClustering/results_cosineSim_weight0.7_decay0.25';

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

%% make main diagonal zero
disp('make main diagonal zero')
high_res_tracts(logical(eye(size(high_res_tracts)))) = 0;

%% load clustering results
disp('load clustering results')
clusters = load(fullfile(results_dir,'clusters.mat'));
high_res_tracts = load(fullfile(results_dir,'high_res_tracts.mat'));
high_res_tracts = high_res_tracts.high_res_tracts;

%% calc SC matrix
disp('calc SC matrix')
numVox = size(clusters.clusterIdPerVoxel, 1);
iterations = [66, 100, 250, 500, 1000];
connmat = cell(length(iterations), 1);
for i=1:length(iterations)
  disp(['i=' num2str(i)])
  iter = iterations(i);
  currentIds = clusters.clusterIdPerVoxel(:,iter);
  
  disp('first loop...')
  connmatTemp = zeros(iter, numVox);
  for clustId = 1:iter
    ids = find( currentIds==clustId );
    connmatTemp(clustId,:) = sum(high_res_tracts(ids,:), 1);
  end
  
  disp('second loop...')
  connmat{i} = zeros(iter,iter);
  for clustId = 1:iter
    ids = find( currentIds==clustId );
    connmat{i}(:,clustId) = sum(connmatTemp(:,ids), 2);
  end
  
end

%% calc ROI size
roi_size = cell(length(iterations), 1);
for i=1:length(iterations)
  roi_size{i} = histc(clusters.clusterIdPerVoxel(:,iterations(i)), 1:iterations(i));
end

%% save
save(fullfile(results_dir,'connmat.mat'), 'iterations', 'connmat', 'roi_size')
