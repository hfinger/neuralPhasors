function   perConn = metricsConn_wd(CIJ)

load('SC4Holger_avg_dist.mat');

% distances between nodes: euclidean, all-pairs-shortest-paths, number steps
perConn.euclDist = avg_dist;
[perConn.apspDist,perConn.stepsDist] = distance_wei(1./CIJ);
corrDist = corr(perConn.euclDist, perConn.apspDist);    

[perConn.EBC, ~] = edge_betweenness_wei(1./CIJ);
