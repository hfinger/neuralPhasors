%% Outlier & Connected Components

% path = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/' RecText '/'];
% 
% 
% OutConn = load([path 'OutlierConnectedComp.mat']);
% Outlier = OutConn.Outliers;
% ConnComp = OutConn.ConnectedComponents;
% 
% 
% % Plotting Outlier
% colMapPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/' RecText '/cmap.mat'];
% cmap = load(colMapPath);
% cmap = cmap.cmap3;
% fig = figure;
% i = imagesc(0:0.1:0.9, 1:1000, cell2mat(Outlier'));
% a = gca;
% colormap(a,cmap);
% title('Sum of Number of Outliers within a cluster over number of clusters and distance weighing');
% xlabel('Degree of Distance Weighing in Composite Matrix');
% ylabel('Number of Clusters');
% c = colorbar;
% ylabel(c, 'Total Number of Outliers');
% makedatatip(i, [0.5 700; 0.5 400; 0.6 550]);
% savefig([path 'NumOutlierWeighingFactor']);
% 
%  java
% % Plotting Connected Components
% colMapPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/cmapmeanconncomp.mat'];
% cmap = load(colMapPath);
% cmap = cmap.cmap;
% fig = figure;
% i = imagesc(0:0.1:0.9, 1:1000, cell2mat(ConnComp'));
% a = gca;
% colormap(a,cmap);
% v = caxis;
% caxis([1 v(2)]);
% title('Mean Number of Connected Components within a cluster over number of clusters and distance weighing');
% xlabel('Degree of Distance Weighing in Composite Matrix');
% ylabel('Number of Clusters');
% c = colorbar;
% ylabel(c, 'Mean Number of Connected Components');
% % makedatatip(i, [0.5 700; 0.5 400; 0.6 550]);
% savefig([path 'MeanConnCompWeighingFactor']);
% 
% %Plotting Connected Components Thresholded
% colMapPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/cmapmeanconncompthresh.mat'];
% cmap = load(colMapPath);
% cmap = cmap.cmap;
% fig = figure;
% i = imagesc(0:0.1:0.9, 1:1000, cell2mat(ConnComp'));
% a = gca;
% colormap(a,cmap);
% caxis([1 1.05]);
% title('Mean Number of Connected Components within a cluster over number of clusters and distance weighing');
% xlabel('Degree of Distance Weighing in Composite Matrix');
% ylabel('Number of Clusters');
% c = colorbar;
% ylabel(c, 'Mean Number of Connected Components');
% % makedatatip(i, [0.5 700; 0.5 400; 0.6 550]);
% savefig([path 'MeanConnCompWeighingFactorThresh']);

%% ClusteringCoefficient


% Clustering Coefficient
statPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1' '/NonRec/StatandMet/'];
statsAll = load([statPath 'statsAll.mat']);
clusCoefMean = statsAll.statsAll.clustering_coefMean;

fig = figure;

clust_coeff = zeros(1000,1000);
for i = 1:1000
    clusCoefMean(i) = clusCoefMean(i)/i;
    if i ==1
        clust_coeff(1:i,i) = NaN;
    else
    clust_coeff(1:i,i) = statsAll.statsAll.clustering_coef{i};
    end
    if i< 1000
    clust_coeff(i+1:1000,i) = NaN;
    end
    clust_coeff(:,i) = clust_coeff(:,i)/i;
    
end
set(0,'DefaultAxesFontSize',18,'DefaultTextFontSize',18);
p = plot(1:1000, clusCoefMean);

q = quantile(clust_coeff, [0.25, 0.5, 0.75], 1);
hold on;
plot(1:1000, q(1,:));
hold on;
plot(1:1000, q(2,:));
hold on;
plot(1:1000, q(3,:))

hold on;
s = nanstd(clust_coeff);

plot(1:1000, s);
    
    title('Mean Clustering Coefficient - Distance Weighting 0', 'fontsize', 15);
    xlabel('Number of Clusters (or Iterations)', 'fontsize', 20);
    ylabel('Clustering Coefficient (log scale)', 'fontsize', 20);
ax = gca;
ax.YScale = 'log';
legend({'Mean Clustering Coefficient', 'First Quantile', 'Median', 'Third Quantile', 'Standard Deviation'}, 'FontSize', 15);
savefig([statPath '/clustcoeffMean0weight.fig']);

%ClusteringCoefficient normalised by number of clusters



