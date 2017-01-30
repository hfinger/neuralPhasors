function clusterStatsNew(clusterPath, clustRange, OutputPath)
Cluster = load(clusterPath);
normClusterConnmat = Cluster.normClusterConnmat;
clusterCoords = Cluster.voxelCoordByCluster;

betwCent = cell(clustRange(end));
Dcell    = cell(clustRange(end));
meanShortestDist1cell = cell(clustRange(end));
meanShortestDist2cell = cell(clustRange(end));
clustering_coefcell = cell(clustRange(end));
lambdacell = cell(clustRange(end));
efficiencycell = cell(clustRange(end));
clusterSizeAll = cell(clustRange(end));
clusterSecondMomentAll = cell(clustRange(end));

% L = cellfun(@(x) 1./x, normClusterConnmat, 'UniformOutput', false);
% % betCentr = cellfun(@(x) betweenness_wei(x), L, 'UniformOutput', false);
% D = cellfun(@(x) distance_wei(x), L, 'UniformOutput', false);
% meanShortestDist1 = cellfun(@(x) mean(x,1), D, 'UniformOutput', false);
% 
% meanShortestDist2 = cellfun(@(x) mean(x,2), D, 'UniformOutput', false);
% [lambda, efficiency] = cellfun(@charpath, D, 'UniformOutput', false);
% clustering_coef = cellfun(@clustering_coef_wu, normClusterConnmat, 'UniformOutput', false);
% clusterSize = cellfun(@(x) length(x), clusterCoords, 'UniformOutput', false);
% clusterSecondMoment = cellfun(@(x) mean(sqrt(sum(bsxfun(@minus,x,mean(x,1)).^2,2))), clusterCoords, 'UniformOutput', false);






for clusterNum = clustRange
    
%    betCentr = nan(clusterNum,1);
%     
    clustConnmat = normClusterConnmat{clusterNum,1};
%     includeClusterIds = find(~isnan(sum(clustConnmat,2)));
%     clustConnmat = clustConnmat(includeClusterIds,includeClusterIds);
%     
    L = 1./(clustConnmat);
%     
%     betCentr(includeClusterIds) = betweenness_wei(L);
%     betwCent{clusterNum} = betCentr;
%     
% 
    D = distance_wei(L);
    Dcell{clusterNum} = D;
    meanShortestDist1 = nan(clusterNum,1);
    meanShortestDist2 = nan(clusterNum,1);
    clustering_coef = nan(clusterNum,1);

%     meanShortestDist1(includeClusterIds) = mean(D,1);
    meanShortestDist1 = mean(D,1);
    meanShortestDist1cell{clusterNum} = meanShortestDist1;
%     meanShortestDist2(includeClusterIds) = mean(D,2);

    meanShortestDist2 = mean(D,2);
    meanShortestDist2cell{clusterNum} = meanShortestDist2;

    [lambda, efficiency] = charpath(D);
    lambdacell{clusterNum} = lambda;
    efficiencycell{clusterNum} = efficiency;
    
%     clustering_coef(includeClusterIds) = clustering_coef_wu(clustConnmat);
    clustering_coef = clustering_coef_wu(clustConnmat);

    clustering_coefcell{clusterNum} = clustering_coef;
    

    clusterSize = cellfun(@(x) length(x), clusterCoords(:,clusterNum), 'UniformOutput', false);
    clusterSecondMoment = cellfun(@(x) mean(sqrt(sum(bsxfun(@minus,x,mean(x,1)).^2,2))), clusterCoords(:,clusterNum), 'UniformOutput', false);
            
    clusterSizeAll{clusterNum} = cell2mat(clusterSize);
    clusterSecondMomentAll{clusterNum} = cell2mat(clusterSecondMoment);
end
    
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end
save ([OutputPath 'clusterstatandmet'], 'betwCent', 'Dcell', ...
    'meanShortestDist1cell', 'meanShortestDistscell', 'clustering_coefcell',...
    'lambdacell', 'efficiencycell', 'clusterSizeAll', 'clusterSecondMomentAll');



end