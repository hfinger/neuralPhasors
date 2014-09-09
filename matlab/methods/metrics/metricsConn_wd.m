function   perConn = metricsConn_wd(CIJ)

load('SC4Holger_avg_dist.mat');
perConn.euclDist = avg_dist;

[perConn.EBC, ~] = edge_betweenness_wei(1./CIJ);
