function [ paths_ordered, pl_ordered ] = getDelayWeightedSWPs( C, threshold, targets, maxPathLength, invC )
%GETDELAYWEIGHTEDSWPS Function that extracts all shortest weighted paths
%between 2 target nodes, weights them according to information transmition
%delays and orderes them from shortest to longest
%
% Input Parameters:
%   C - 2d connectivity matrix (will be binarized with threshold)
%   threshold - Connectivity threshold used to binarize C (scalar)
%   targets - vector with 2 indices indicating the two nodes to search for paths betwenn
%   maxPathLength - maximum number of edges for each path between the 2 nodes
%   invC - if True, pathlengths will be calculated based on 1./C (should be
%          used for weighted Cs)
%
% Output:
% paths_ordered - cell with one entry for each SWP (ordered in increasing length
% pl_ordered - summed up pathlengths for each path in paths_ordered

%% find all paths of a certain maximum lenght between target nodes, get their connection weights and time delays and create delay-weighted shortest-weighted paths

% get all paths from target node 1 of specified maximum length
[~,~,~,~,allpaths,~] = findpaths(C > threshold, targets(1), maxPathLength, 1);

% extract paths that include specified target node
n = 1;
for p=1:size(allpaths,2)
    p_tmp = allpaths(:,p);
    if sum(p_tmp == targets(2)) > 0
        idx = find(p_tmp == targets(2));
        paths{n} = p_tmp(1:idx);
        n = n+1;
    end
end

% get rid of path multiples
finalPaths = {};
n = 1;
for p=1:length(paths)
    inc = true;
    for p2 = 1:length(finalPaths)
        if length(paths{1,p}) == length(finalPaths{1,p2})
            if paths{1,p} == finalPaths{1,p2}
                inc = false;
            end
        end
    end
    if inc
        finalPaths{n} = paths{1,p};
        n = n+1;
    end
end
paths = finalPaths;

% calculate pathlengths based on connectivity matrix
if invC
    C = 1./C;
end
pathlengths = zeros(1,length(paths));
for p=1:length(paths)
    for i=1:length(paths{1,p})-1
        pathlengths(p) = pathlengths(p) + C(paths{1,p}(i),paths{1,p}(i+1));
    end
end

% sort pathlengths
[pl_ordered,idx] = sort(pathlengths);

% order paths according to indices
for i=1:length(idx)
    paths_ordered{i} = paths{1,idx(i)};
end

end

