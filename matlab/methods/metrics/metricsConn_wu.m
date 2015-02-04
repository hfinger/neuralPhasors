function   perConn = metricsConn_wu(CIJ)

paths = dataPaths();
load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));

perConn.euclDist = avg_dist;

[perConn.EBC, ~] = edge_betweenness_wei(1./CIJ);