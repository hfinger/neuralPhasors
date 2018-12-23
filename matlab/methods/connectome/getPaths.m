function [node_pairs,shortest_path_length,paths_all] = getPaths(C,pairs1,pairs2,maxPathLength,minPathNum)
%GETPATHS Summary of this function goes here
%   Detailed explanation goes here

% find paths between each pair of nodes (non-symmetric). store paths, 
% shortest path lengths and pairs
nodes = 1:size(C,1);
paths_all = cell(0);
node_pairs = cell(0);
shortest_path_length = cell(0);
n = 1;
for i=1:length(nodes)
    [~,~,~,~,paths,~] = findpaths(C > 0, i, maxPathLength, 1); % find all paths
    for j=1:length(nodes)
        pathsTmp = paths;
        idx_col = find(sum(pathsTmp == j,1) > 0); % find all paths from node i to node j
        
        if j == i || isempty(idx_col) % disregard self connections
            continue
        end
        
        if minPathNum > 1
            if C(i,j) > 0 || C(j,i) > 0 % disregard direct
                continue
            end
        end
        
        idx_row = zeros(1,length(idx_col));
        for k=1:length(idx_col) 
            idx_row(k) = find(pathsTmp(:,idx_col(k)) == j); % find position of node j on path
            pathsTmp(idx_row(k)+1:end,idx_col(k)) = 0; % set all following nodes to 0
        end
        
        debugTest = sum(pathsTmp(:,idx_col)==i,2);
        if sum(debugTest(2:end)) > 0
            disp('error: source node is within path!!')
        end
        
        paths_all{1,n} = pathsTmp(:,idx_col);
        shortest_path_length{n} = min(idx_row) - 1;
        node_pairs{n,1} = [i,j];
        n = n+1;
    end
end

node_pairs = cell2mat(node_pairs);
for i=1:size(node_pairs,1)
    paths_all{i} = unique(paths_all{i}','rows')';
end

% check for bi-directional double entries
% for i=1:size(node_pairs,1)
%     paths = unique(paths_all{i}','rows')';
%     if size(paths,2) < minPathNum
%         node_pairs(i,:) = 0; % disregard all node pairs with less than minPathNum connections
%     end
%     idx = all(node_pairs == repmat(fliplr(node_pairs(i,:)),size(node_pairs,1),1),2);
%     node_pairs(idx,:) = 0; % disregard flipped node pairs
% end


idx = node_pairs(:,1) > 0;
node_pairs = node_pairs(idx,:);
paths_all = {paths_all{idx}};
shortest_path_length = {shortest_path_length{idx}};
shortest_path_length = cell2mat(shortest_path_length);

end

