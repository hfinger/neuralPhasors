close all;
clear all;
clc;

savepath_prefix = '/net/store/nbp/projects/phasesim/workdir/Holger/20160330_GraclusClustering/';

%% load tractography data:
disp('load tractography data')
path_prefix = '/net/store/nbp/projects/phasesim/databases/DTI_subject01/';

high_res_tracts = load([path_prefix 'tractographyData/high_res_tracts.mat']);
high_res_tracts = high_res_tracts.high_res_tracts;

distance_matrix = load('/net/store/nbp/projects/phasesim/src_enes2/code&data/distances/250_distance_matrix');
distance_matrix = distance_matrix.distance_matrix;

%% make main diagonal zero
disp('make main diagonal zero')
high_res_tracts = high_res_tracts - diag(diag(high_res_tracts));

%% convert to probabilities:
disp('convert to probabilities')
high_res_tracts = bsxfun(@rdivide,high_res_tracts,sum(high_res_tracts,2));

%% make the matrix symmetrical
disp('make the matrix symmetrical')
high_res_tracts = high_res_tracts + high_res_tracts';

%% convert to cosine similarity matrix:
useCosineSimilarity = false;
if useCosineSimilarity
  disp('convert to cosine similarity matrix')
  tmp = high_res_tracts*high_res_tracts';
  tmp = bsxfun(@rdivide,tmp, sqrt(sum(high_res_tracts.^2,2)) );
  high_res_tracts = bsxfun(@rdivide,tmp, sqrt(sum(high_res_tracts.^2,1)) );
  clear tmp;
  
  % reset main diagonal to zero:
  disp('reset main diagonal to zero')
  high_res_tracts = high_res_tracts - diag(diag(high_res_tracts));
end

%% calculate exponential decay for the distance matrix
disp('calculate exponential decay for the distance matrix')
decay_constant = 1;
distance_matrix = spfun(@exp, -distance_matrix * decay_constant);

%% normalize the matrices:
disp('normalize the matrices')
% distance_matrix = bsxfun(@rdivide, distance_matrix, sum(distance_matrix,2));
% high_res_tracts = bsxfun(@rdivide, high_res_tracts, sum(high_res_tracts,2));
% distance_matrix = bsxfun(@rdivide, distance_matrix, mean(distance_matrix(:)));
% high_res_tracts = bsxfun(@rdivide, high_res_tracts, mean(high_res_tracts(:)));
distance_matrix = size(high_res_tracts,1) * bsxfun(@rdivide, distance_matrix, sum(distance_matrix(:)));
high_res_tracts = size(high_res_tracts,1) * bsxfun(@rdivide, high_res_tracts, sum(high_res_tracts(:)));

%% combine DTI data with distance matrix:
disp('combine DTI data with distance matrix')
distance_weight = 0.5;
connectivity_weight = 1 - distance_weight;
high_res_tracts = (distance_weight*distance_matrix) + (connectivity_weight*high_res_tracts);

%% sparsity stats:
totalNumElem = length(high_res_tracts(:));
numNonZero = sum(high_res_tracts(:)>0);
numAboveThresh = sum(high_res_tracts(:)>0.0005);
numBelowThreshAndNonZero = numNonZero-numAboveThresh;
numBelowThreshAll = totalNumElem-numAboveThresh;
save('sparsityStatsSumNormTimesNumRows.mat','totalNumElem','numNonZero','numAboveThresh','numBelowThreshAndNonZero','numBelowThreshAll')

%% do clustering
disp('do clustering')
recursive_splitting = false;
cluster_count = 500;
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
results_dir = ['new_results_' results_dir];
results_path = fullfile(savepath_prefix, results_dir);
mkdir(results_path)
save(fullfile(results_path,'clusters.mat'), 'clusterIdPerVoxel', 'largestClusterId', 'cutValue')
save(fullfile(results_path,'high_res_tracts.mat'), 'high_res_tracts', '-v7.3')


