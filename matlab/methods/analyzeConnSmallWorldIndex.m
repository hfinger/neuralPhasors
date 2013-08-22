function [ stats ] = analyzeConnSmallWorldIndex( W )
%ANALYZECONNSMALLWORLDINDEX Summary of this function goes here
%   Detailed explanation goes here

fids=sort(unique(W.f0));

% resample connections but with same radius:
connRadius = sqrt(W.dx.^2+W.dy.^2);
connAngle = 2*pi*rand(size(connRadius));
WRandom.dx = round(connRadius.*cos(connAngle));
WRandom.dy = round(connRadius.*sin(connAngle));
WRandom.f1 = randi(max(fids),size(W.f1));
WRandom.f0 = W.f0;
WRandom.w = W.w;

disp('calc shortest-path of real network')
stats.ShortestPath = analyzeConnShortestPath( W, 50 );
disp('calc cluster-index of real network')
stats.ClusterIndex = analyzeConnClusterIndex( W );

% compute again for random connectivity:
disp('calc shortest-path of random network')
stats.ShortestPathRandom = analyzeConnShortestPath( WRandom, 50 );
disp('calc cluster-index of random network')
stats.ClusterIndexRandom = analyzeConnClusterIndex( WRandom );

% calculate global statistics:
stats.meanShortestPath = mean(stats.ShortestPath.shortestPathLengths(:));
stats.meanClusterIndex = mean(stats.ClusterIndex.clustering_coeffs(:));
stats.meanShortestPathRandom = mean(stats.ShortestPathRandom.shortestPathLengths(:));
stats.meanClusterIndexRandom = mean(stats.ClusterIndexRandom.clustering_coeffs(:));

% calculate small world index:
stats.SmallWorldIndex = (stats.meanClusterIndex / stats.meanClusterIndexRandom) / (stats.meanShortestPath / stats.meanShortestPathRandom);

stats.WRandom = WRandom;

end

