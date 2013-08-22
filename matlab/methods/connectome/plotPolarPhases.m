function plotPolarPhases( localPhase, numCluster, N )
%PLOTPOLARPHASES Summary of this function goes here
%   Detailed explanation goes here

cmap=colormap('lines');
% clf;
for clust=1:numCluster
  ids=(clust-1)*N+1:clust*N;
  h = polar(repmat(localPhase(ids),1,2)',[zeros(size(localPhase(ids))) ones(size(localPhase(ids)))]');
%   hold on;
  hold all; 
  set(h,'Color',cmap(clust,:))
end

end

