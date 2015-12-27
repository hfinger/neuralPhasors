function   perConn = metricsConn_wu(CIJ)

paths = dataPaths();
load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));

avg_dist(isnan(avg_dist)) = 0;
perConn.euclDist = avg_dist;                                                % euclidean distance between all pairs of ROI 

[perConn.floydDist, perConn.floydNoEdges] = distance_wei(1./CIJ);           % all-pairs shortest path distances and
                                                                            % numbers of edges per shortest path
                                                                            
                                                                            % notice that: mean(perConn.floydDist)
                                                                            % equals the characteristic path length,
                                                                            % i.e. perGlob.lambda in metricsGlobal_wu

[perConn.EBC, ~] = edge_betweenness_wei(1./CIJ);                            % weighted edge betweenness centrality:
                                                                            % percentage of shortest paths containing a given edge