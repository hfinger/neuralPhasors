%%
Y=simResult.Y;
numCluster=network.numCluster;
N=network.N;

%%
globPhase = angle(sum(exp(1i*Y),1));
localPhase = mod(bsxfun(@minus,Y,globPhase),2*pi);
clusterPhase = zeros(numCluster,size(localPhase,2));
clusterOrderCoeff = zeros(numCluster,size(localPhase,2));
for clust=1:numCluster
  ids=(clust-1)*N+1:clust*N;
  clusterPhase(clust,:) = angle(sum(exp(1i*Y(ids,:)),1)); 
  clusterOrderCoeff(clust,:) = abs(sum(exp(1i*Y(ids,:)),1)/size(Y,1)); 
end

%%
figure(1); clf;
imagesc(sin(Y))

figure(2); clf;
imagesc(1000*diff(Y,[],2)/(2*pi));
xlabel('time [ms]')

colormap jet;
colorbar;

% figure(3); clf;
% imagesc(localPhase); 
% colormap hsv; 
% colorbar

figure(4); clf;
imagesc(C);
xlabel('source neuron')
ylabel('target neuron')

% figure(5); clf;
% imagesc(clusterPhase)

% figure(6); clf;
% plotPolarPhases( localPhase(:,1), numCluster, N );

figure(8);
plot(clusterOrderCoeff')
xlabel('time [ms]')
clusterNames=num2cell(num2str((1:size(clusterOrderCoeff,1))'));
legend(clusterNames{:})