 datapaths = dataPaths();       
for subjNum = 1:22
            clusterPath = [datapaths.workdir '/Arushi/20150824Allsubjectrecursivencut/'];
            outputPath = [datapaths.workdir '/Arushi/20150920StatsAndMetrics/'];
               clusterType = 'fullcos';
                
     clusterNum = 1001;
                    %             clusterNum = mod(this.params.ProbClustbetcent.split, 1000) + 1;
                    %             if clusterNum == 1
                    %                 clusterNum = (clusterNum * 1000) + 1;
                    %             end
                    
                  
                        if exist([clusterPath num2str(subjNum)], 'dir')
                            disp(['subj' num2str(subjNum) clusterType num2str(clusterNum)]);

                            subjClusterPath = [clusterPath num2str(subjNum) '/' clusterType '/' clusterType '/postprocessing/'];
                            
                            disp(['subj' num2str(subjNum) 'cluster' num2str(clusterNum)]);
                            filepath = [subjClusterPath 'clusterConnmat' num2str(clusterNum) '.mat'];
                            clusterConnmat = load(filepath);
                            clusterConnmat = clusterConnmat.clusterConnmat;
                            clusterConnmat = bsxfun(@rdivide, clusterConnmat, sum(clusterConnmat,2));
                            
                            includeClusterIds = find(~isnan(sum(clusterConnmat,2)));
                            
                            connSym = clusterConnmat + clusterConnmat';
                            connSym = connSym(includeClusterIds,includeClusterIds);
                            connSym(logical(eye(size(connSym)))) = 0;
                            L = 1./(connSym);
                            
                            clustering_coef = nan(size(clusterConnmat,1),1);
                            meanShortestDist1 = nan(size(clusterConnmat,1),1);
                            meanShortestDist2 = nan(size(clusterConnmat,1),1);

                            clustering_coef(includeClusterIds) = clustering_coef_wu(connSym);
                            D = distance_wei(L);
                            meanShortestDist1(includeClusterIds) = mean(D,1);
                            meanShortestDist2(includeClusterIds) = mean(D,2);
                            [lambda, efficiency] = charpath(D);
           
                            newOutputPath = [outputPath num2str(subjNum) '/' clusterType '/'];
                            if ~exist(newOutputPath, 'dir')
                                mkdir(newOutputPath);
                            end
           
                            save([newOutputPath clusterType num2str(clusterNum) '.mat'], 'meanShortestDist1', 'meanShortestDist2', 'lambda', 'efficiency', 'clustering_coef');
                                  
                        end
                    end
                    %%%% END EDIT HERE:                                %%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            
       