% Generate connectivity matrices from tracking data for all subjects, cluster them and analyse them
% 20160407 Added CompSimMatCalc

%% Parameters
fn = 18;
%% Get Full Connectivity and Waytotal matrices from grid job results of ProbtrackX
if fn == 1
    %Parameters
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160415_Probtrackxallsubjects';
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160417_Allconnmat'];
    subjRange = [6];%[1:4, 6:13, 15, 17:22];
    splitPerSubject = 100;
    
    %Call Function
    GetFullConnmat(inputPath, outputPath, subjRange, splitPerSubject);
end

%% Get FS Connectivity and Waytotal Matrices from full connectivity matrix above
if fn == 2
    %Parameters
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150821_FullConnmat';
    inputMaskPath = '/net/store/nbp/projects/phasesim/databases/Bastian_DTI/masks/';
    inputcomplFSPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_FSConnmat/'];
    subjTotal = 22;
    
    %Call Function
    GetFSConnmat(inputPath, inputMaskPath, inputcomplFSPath, outputPath, subjTotal);
end
%% Remove overlapping FS voxels for all subjects
if fn == 3
    %Parameters
    inputMaskPath = '/net/store/nbp/projects/phasesim/databases/Bastian_DTI/masks/';
    inputcomplFSPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_overlappingFS/'];
    subjTotal = 22;
    
    %Call Function
    removeOverlappingFS(inputMaskPath, inputcomplFSPath, outputPath, subjTotal);
    
end
%% Get Recursive Ncut from full connectivity matrix
if fn == 4
    
    %Parameters
    useCosineSimilarity = 0;
    useFSROI = 1;
    numberIterations = 1000;
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150821_FullConnmat';
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/'];
    clusterPath =   '/net/store/nbp/projects/phasesim/workdir/Arushi/20150824_overlappingFS/';
    subjTotal = 22;
    
    
    %Call Function
    GetRecursiveNcut(inputPath, clusterPath, useCosineSimilarity, useFSROI, numberIterations, outputPath, subjTotal);
    
end

%% Ncut post processing

if fn == 5
    subjTotal = 22;
    clusterTotal = 1000;
    clusterPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/';
    clusterType = 'fullconn';
    maskPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
    connmatPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150821_FullConnmat';
    
    %CallFunction
    RecursiveNcutPostProcess(subjTotal, clusterTotal, clusterPath, clusterType, maskPath, connmatPath);
end

%% Statistics and Metrics

if fn==6
    
    subjTotal = 22;
    clusterTotal = 1000;
    clusterPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/';
    statAndMetricPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/201509100Allsubjectstats';
    clusterType = 'fullcos';
    
    %CallFunction
    statsAndMetrics(subjTotal, clusterTotal, clusterType, clusterPath, statAndMetricPath);
end
%** Note - It was ultimately run in grid job. Find program in
%matlab/classes/ProbClustmetrics for other metrics
%For betweenness centrality - matlab/classes/ProbClustbetcent
%% Cluster Stats

if fn == 7
    
    %CallFunction
    clusterStats();
end

%% Collect subject stats and metrics to one file

if fn == 8
    filepath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150920StatsAndMetrics/';
    subjRange = 1:11;
    
    collectstatsandmetrics(filepath, subjRange);
end

%% Plot subjects stats and metrics

if fn == 9
    
    filepath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150920StatsAndMetrics/';
    subjRange = 1:11;
    
    plotstatsandmetrics(filepath, subjRange);
end

%% Make clusternii

if fn == 10
    
    clusterfilepath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/';
    maskfilepath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
    savepath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20151005clustermask/';
    subjRange = 1;
    modeRange = 3;
    clustRange = 2:5;
    
    makeclusternii(clusterfilepath, maskfilepath, savepath, modeRange, subjRange, clustRange);
end

%% Composite Similarity Matrix
% 20160407

if fn == 11
    
    run(CompSimMatCalcFn);
    
end

%% Graclus recursive Ncut

if fn == 12
    
    run(graclusNcut);
    
end

%% Global cut value calculation
if fn == 13
    run(GlobalCutValCalcscript);
    
end

%% Collect Graclus Result for NonRecursive Clusters into one matrix\

if fn ==14
    
    run(collectGraclusNonRecResult);
end

%% Generate clusterwise connmat for all number of clusters
if fn ==15
    run(GenerateClusterConnmatScript);
    
end
%% Calculate cluster stats


if fn == 16
    run(clusterStatsNewScript);
end


%% Calculate intersecting clusters

if fn == 17
    run(calcCommonClustersScript);
    subjRange = 1;
    clustRange = 2:1000;
    for subjNum = subjRange
        for clusterNum = clustRange
            OutputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/NonRec/IntersectingClusters';
            voxelIndByCluster1 = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/NonRec/Clustorg/detailsClust'...
                num2str(clusterNum)]);
            voxelIndByCluster1 = voxelIndByCluster1.voxelIndByCluster;
            voxelIndByCluster2 = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/NonRec/Clustorg/detailsClust'...
                num2str(clusterNum+1)]);
            voxelIndByCluster2 = voxelIndByCluster2.voxelIndByCluster;
            
            calcCommonClusters(voxelIndByCluster1, voxelIndByCluster, OutputPath, clusterNum);
        end
    end
end

%% Calculate outliers

if fn == 18
    clear all;
    subjRange = 1;
    clustRange = 2:1000;
    RecursiveText = 'Rec';
    cosText = 'conn';
    WeighingFacRange = 0:0.1:1;
    for subjNum = subjRange
        for WeighingFactor = WeighingFacRange
            load('/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/FinalCoord/FinalCoord1.mat');
            clustPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_GraclusCut/',...
                cosText '/' RecursiveText '/decay-1weigh' num2str(WeighingFactor) '/1/normbysumthresh1/'];
            switch RecursiveText
                case 'Rec'
                    load([clustPath 'graclusResultClust1000.mat']);
                    number_of_outliers = zeros(1000,1000);
                    number_of_connectedComponents = zeros(1000,1000);
                    number_of_outliers(:,1) = NaN;
                    number_of_connectedComponents(:,1) = NaN;
                    for clusterNum = clustRange
                        a = stepWiseClusters(:,clusterNum);
                        number_of_outliers(clusterNum+1:end,clusterNum) = NaN;
                        number_of_connectedComponents(clusterNum+1:end,clusterNum) = NaN;
                        
                        
                        %                 detClust = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW'...
                        %                     num2str(WeighingFactor) '/detailsClust' num2str(clusterNum) '.mat']);
                        
                        for singleCluster = 1:clusterNum
                            %                     voxel_coords = voxelCoordByCluster{singleCluster};
                            voxel_coords = FinalCoord((a==singleCluster), :);
                            [number_of_outliers(singleCluster, clusterNum), number_of_connectedComponents(singleCluster, clusterNum),clusterSecondMoment] = connectedComp(voxel_coords);
                        end
                    end
                    OutputPath = [PostProcessPath RecursiveText '/decay' num2str(decayParam)...
                        'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh'...
                        num2str(threshFactor) '/' cosText '/' 'Subj'...
                        num2str(subjNum) '/'];
                    save([OutputPath 'components.mat'], 'number_of_outliers', 'number_of_components', '-v7.3');
                    
                case 'NonRec'
                    number_of_outliers = zeros(1000,1000);
                    number_of_connectedComponents = zeros(1000,1000);
                    number_of_outliers(:,1) = NaN;
                    number_of_connectedComponents(:,1) = NaN;
                    for clusterNum = clustRange
                        load([clustPath 'graclusResultClust' num2str(clusterNum) '.mat']);
                        a = stepWiseClusters;
                        number_of_outliers(clusterNum+1:end,clusterNum) = NaN;
                        number_of_connectedComponents(clusterNum+1:end,clusterNum) = NaN;
                       
                        for singleCluster = 1:clusterNum
                             voxel_coords = FinalCoord((a==singleCluster), :);
                        [number_of_outliers(singleCluster, clusterNum), number_of_connectedComponents(singleCluster, clusterNum),clusterSecondMoment] = connectedComp(voxel_coords);
                        
                        end
                        
                        
                    end
                     OutputPath = [PostProcessPath RecursiveText '/decay' num2str(decayParam)...
                        'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh'...
                        num2str(threshFactor) '/' cosText '/' 'Subj'...
                        num2str(subjNum) '/'];
                    save([OutputPath 'components.mat'], 'number_of_outliers', 'number_of_components', '-v7.3');
            end
        end
    end
end
%% Reduce stats and metrics to one file
    if fn == 19
        clusterRange = 2:1000;
        
        for clusterNum = clusterRange
            StatTemp= load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW0/StatandMet/clusterstatandmetClust' num2str(clusterNum) '.mat']);
            statsAll.numROI{clusterNum} = clusterNum;
            statsAll.meanShortestDist1{clusterNum} = StatTemp.meanShortestDist1;
            statsAll.meanShortestDist2{clusterNum} = StatTemp.meanShortestDist2;
            statsAll.lambda{clusterNum} = StatTemp.lambda;
            statsAll.efficiency{clusterNum} = StatTemp.efficiency;
            statsAll.clustering_coef{clusterNum} = StatTemp.clustering_coef;
            statsAll.clusterSecondMoment{clusterNum} = cell2mat(StatTemp.clusterSecondMoment);
            statsAll.clusterSize{clusterNum} = cell2mat(StatTemp.clusterSize);
            statsAll.D{clusterNum} = StatTemp.D;
            
        end
        statsAll.numROI = statsAll.numROI';
        statsAll.meanShortestDist1 = statsAll.meanShortestDist1';
        statsAll.meanShortestDist2 = statsAll.meanShortestDist2';
        statsAll.lambda = statsAll.lambda';
        statsAll.efficiency = statsAll.efficiency';
        statsAll.clustering_coef = statsAll.clustering_coef';
        statsAll.clusterSecondMoment = statsAll.clusterSecondMoment';
        statsAll.clusterSize = statsAll.clusterSize';
        statsAll.D =  statsAll.D';
        statsAll.clusterSecondMomentMean = cellfun(@nanmean, statsAll.clusterSecondMoment);
        
        statsAll.clusterSecondMomentStd = cellfun(@nanstd, statsAll.clusterSecondMoment);
        statsAll.clusterSizeMean = cellfun(@nanmean, statsAll.clusterSize);
        statsAll.clustering_coefMean = cellfun(@nanmean, statsAll.clustering_coef);
        statsAll.efficiency = cell2mat(statsAll.efficiency);
        statsAll.lambda = cell2mat(statsAll.lambda);
        statsAll.numROI = cell2mat(statsAll.numROI);
        
        statsAll.numROI(2:1000) = statsAll.numROI(1:999);
        statsAll.numROI(1) = 1;
        
        statsAll.efficiency(2:1000) = statsAll.efficiency(1:999);
        statsAll.efficiency(1) = NaN;
        
        statsAll.lambda(2:1000) = statsAll.lambda(1:999);
        statsAll.lambda(1) = NaN;
        save('/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW0/StatandMet/statsAll', 'statsAll', '-v7.3');
    end
    
    %% Plot Metrics
    if fn == 20
        
        yname = 'clustering_coefMean';
        figure(1)
        clf;
        hold on;
        
        plot(statsAll.numROI, statsAll.(yname))
        
        ylabel('clustering coef')
        xlabel('resolution [#ROI]')
        % set(gca,'ylim',[1e-4 1e-2])
        % set(gca,'YScale','log');
        saveas(gcf, ['figures/mean_' yname],'png');
    end
    
    %% Processing Outlier and connected components data
    if fn == 21
        ConnectedComponents = cell(10,1);
        Outliers = cell(10,1);
        for WeighingFactor = 0:0.1:0.9
            NumberConnectedComponentsWithin = zeros(1000,1);
            NumberOutliers = zeros(1000,1);
            for clusterNum = 2:1000
                Path = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW' num2str(WeighingFactor) '/Outliers/outlier' num2str(clusterNum)];
                outlier = load(Path);
                NumberConnectedComponentsWithin(clusterNum) = mean(outlier.number_of_connectedComponents);
                NumberOutliers(clusterNum) = sum(outlier.number_of_outliers);
            end
            ConnectedComponents{(round(WeighingFactor*10) +1)} = NumberConnectedComponentsWithin;
            Outliers{(round(WeighingFactor*10)+1)} = NumberOutliers;
        end
        
        save('/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/OutlierConnectedComp', 'ConnectedComponents', 'Outliers');
    end
    
    %% Plotting statistics
    if fn == 22
        run(plotStat);
    end
    %% Construct VoxelByInd cell structure for FS masks
    if fn == 23
        run(makeVoxelByIndForFS);
    end
    
    %% Betweenness Centrality
    
    if fn == 24
        run(ProbClustbetcentNewScript);
    end
    
    %% Combine CalcIntersecting Data
    
    if fn == 25
        
        run(combineCalcIntersectingScript);
    end

    
    %% Combine CalcIntersecting Data
    
    if fn == 26
        
        run(createClusteringVideo);
    end
