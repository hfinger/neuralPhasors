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
p = 2;
v = 4.5;
C = C + 0.1 * diag(ones(size(C,1)/2,1),size(C,1)/2) + 0.1 * diag(ones(size(C,1)/2,1),-size(C,1)/2);
C = bsxfun(@rdivide,C,sum(C.^p,2).^(1/p));
D = D/(v*1e03);

%% find all paths of a certain maximum lenght between target nodes, get their connection weights and time delays and create delay-weighted shortest-weighted paths

% get all paths from target node 1 of specified maximum length
stimPos = [4,31];
maxPathLength = 4;
[~,~,~,~,allpaths,~] = findpaths(C > 0.05, stimPos(1), maxPathLength, 1);

% extract paths that include specified target node
n = 1;
for p=1:size(allpaths,2)
    p_tmp = allpaths(:,p);
    p_tmp = p_tmp(p_tmp ~= 0);
    if sum(p_tmp == 31) > 0
        idx = find(p_tmp == 31);
        paths{n} = p_tmp(1:idx);
        n = n+1;
    end
end

% calculate pathlengths based on delay-weighted connectivity matrix
C_dw = C ./ D;

pathlengths = zeros(1,length(paths));
for p=1:length(paths)
    for i=1:length(paths{1,p})-1
        pathlengths(p) = pathlengths(p) + C(paths{1,p}(i),paths{1,p}(i+1));
    end
end

% sort pathlengths
[pl_ordered,idx] = sort(pathlengths);
pl_ordered = fliplr(pl_ordered);
idx = fliplr(idx);

% order paths according to indices
for i=1:length(idx)
    paths_ordered{i} = paths{1,idx(i)};
end

%% plot k shortest ordered paths in brain

k = 20;

% set node and edge values
edgeVals = zeros(size(C));
for i=1:k
    for j=1:length(paths_ordered{1,i})-1
        edgeVals(paths_ordered{1,i}(j),paths_ordered{1,i}(j+1)) = C_dw(paths_ordered{1,i}(j),paths_ordered{1,i}(j+1));
    end
end
nodeVals = zeros(1,size(C,1));
nodeVals(stimPos) = 1;

% plot nodes and edges in brain
edgeMin = 0.001;
edgeVals = edgeVals * 0.05;
plotBrainConnectivity(nodeVals,edgeVals,edgeMin,1,0)
