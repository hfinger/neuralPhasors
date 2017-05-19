function clusterStatsNew(clustRange,clusterPath, OutputPath)
            clustering_coef = nan(clustRange(end));
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
           
                Cluster = load([clusterPath 'detailsClust2to1000.mat']);
                normClusterConnmat = Cluster.normClusterConnmat;
                clusterCoords = Cluster.voxelCoordByCluster;
                clusterSizeAll = cellfun(@(x) length(x), clusterCoords, 'UniformOutput', false);
                clusterSecondMomentAll = cellfun(@(x) mean(sqrt(sum(bsxfun(@minus,x,mean(x,1)).^2,2))), clusterCoords, 'UniformOutput', false);

            for splitNum = clustRange
                 disp(num2str(splitNum));
                
                %    betCentr = nan(clusterNum,1);
                %
                clustConnmat = normClusterConnmat{splitNum};
                %     includeClusterIds = find(~isnan(sum(clustConnmat,2)));
                %     clustConnmat = clustConnmat(includeClusterIds,includeClusterIds);
                %
                
                
                %     clustering_coef(includeClusterIds) = clustering_coef_wu(clustConnmat);
                
                clustering_coef(1:splitNum, splitNum) = clustering_coef_wu(clustConnmat);
                
                
%                 clusterSize = cellfun(@(x) length(x), clusterCoords(1:split,split), 'UniformOutput', false);
                
%                {split} = cell2mat(clusterSize);
            end
            
            if ~exist(OutputPath, 'dir')
                mkdir(OutputPath);
            end
            save ([OutputPath 'clusterstatandmet'],...
                'clustering_coef', 'clusterSizeAll', 'clusterSecondMomentAll', '-v7.3');