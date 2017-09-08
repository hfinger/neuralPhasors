function [ paths_ordered, pl_ordered ] = getDelayWeightedSWPs( targets, maxPathLength, p, v  )
%GETDELAYWEIGHTEDSWPS Function that extracts all shortest weighted paths
%between 2 target nodes, weights them according to information transmition
%delays and orderes them from shortest to longest
%
% Input Parameters:
%   targets - vector with 2 indices indicating the two nodes to search for paths betwenn
%   maxPathLength - maximum number of edges for each path between the 2 nodes
%   p - order of the norm with which to normalize the connectivity matrix
%   v - velocity of information transfer in the network [m/s]
%
% Output:
% paths_ordered - cell with one entry for each SWP (ordered in increasing length
% pl_ordered - summed up pathlengths for each path in paths_ordered

%% Load connectivity/distance matrix, resort areas and normalize them

% loading
path_SCmat = '/net/store/nbp/projects/phasesim/databases/avg_SC.mat';
path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
load(path_SCmat);
load(path_ResortIDs);
resortIds = [resortIdsMarlene, resortIdsMarlene + 33];

% resorting
C = C(resortIds,:);
C = C(:,resortIds);
D = D(resortIds,:);
D = D(:,resortIds);

%normalizing
C = C + 0.1 * diag(ones(size(C,1)/2,1),size(C,1)/2) + 0.1 * diag(ones(size(C,1)/2,1),-size(C,1)/2);
C = bsxfun(@rdivide,C,sum(C.^p,2).^(1/p));
D = D/(v*1e03);

%% find all paths of a certain maximum lenght between target nodes, get their connection weights and time delays and create delay-weighted shortest-weighted paths

% get all paths from target node 1 of specified maximum length
[~,~,~,~,allpaths,~] = findpaths(C > 0.05, targets(1), maxPathLength, 1);

% extract paths that include specified target node
n = 1;
for p=1:size(allpaths,2)
    p_tmp = allpaths(:,p);
    p_tmp = p_tmp(p_tmp ~= 0);
    if sum(p_tmp == targets(2)) > 0
        idx = find(p_tmp == targets(2));
        paths{n} = p_tmp(1:idx);
        n = n+1;
    end
end

% calculate pathlengths based on delay-weighted connectivity matrix
C_dw = C ./ D;

pathlengths = zeros(1,length(paths));
for p=1:length(paths)
    for i=1:length(paths{1,p})-1
        pathlengths(p) = pathlengths(p) + C_dw(paths{1,p}(i),paths{1,p}(i+1));
    end
    pathlengths(p) = pathlengths(p) / (length(paths{1,p})-1);
end

% sort pathlengths
[pl_ordered,idx] = sort(pathlengths);
pl_ordered = fliplr(pl_ordered);
idx = fliplr(idx);

% order paths according to indices
for i=1:length(idx)
    paths_ordered{i} = paths{1,idx(i)};
end


end

