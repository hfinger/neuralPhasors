function statsAndMetrics( subjTotal, clusterTotal, clusterType, clusterPath, statAndMetricPath )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    if ~subjTotal
        subjTotal = 22;
    end
    
    if ~clusterTotal
        clusterTotal = 1000;
    end
    
    if ~clusterPath
        clusterPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/';
    end
    
    if ~statAndMetricPath
        statAndMetricPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/201509100Allsubjectstats/';
    end
    
    disp(clusterType);
    
    if clusterType == 'fsconn'|| 'fscos'
        startCluster = 67;
    else
        startCluster = 2;
    end

   for subjNum = 1:subjTotal
       if ~exist([clusterPath num2str(subjNum)], 'dir')
           continue
       end
       
       subjClusterPath = [clusterPath num2str(subjNum) '/' clusterType '/' clusterType '/postprocessing/'];
       
       for clusterNum = startCluster:clusterTotal
           filepath = [subjClusterPath 'clusterConnmat' num2str(clusterNum) '.mat'];
           clusterConnmat = load(filepath);
           clusterConnmat = clusterConnmat.clusterConnmat;
           clusterConnmat = bsxfun(@rdivide, clusterConnmat, sum(clusterConnmat,2));
           
           includeClusterIds = find(~isnan(sum(clusterConnmat,2)));

           connSym = clusterConnmat + clusterConnmat';
           connSym = connSym(includeClusterIds,includeClusterIds);
           connSym(logical(eye(size(connSym)))) = 0;
           
           clustering_coef = nan(size(clusterConnmat,1),1);
           betCentr = nan(size(clusterConnmat,1),1);
           meanShortestDist1 = nan(size(clusterConnmat,1),1);
           meanShortestDist2 = nan(size(clusterConnmat,1),1);

           
           clustering_coef(includeClusterIds) = clustering_coef_wu(connSym);
           
           betCentr(includeClusterIds) = betweenness_wei(L);
           
           D = distance_wei(L);
           meanShortestDist1(includeClusterIds) = mean(D,1);
           meanShortestDist2(includeClusterIds) = mean(D,2);
           [lambda, efficiency] = charpath(D);
       end
   end
end

