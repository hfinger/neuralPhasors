% Generate connectivity matrices from tracking data for all subjects, cluster them and analyse them
% 20160407 Added CompSimMatCalc

%% Parameters
run = 18;
%% Get Full Connectivity and Waytotal matrices from grid job results of ProbtrackX
if run == 1
    %Parameters
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160415_Probtrackxallsubjects';
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160417_Allconnmat'];
    subjRange = [6];%[1:4, 6:13, 15, 17:22];
    splitPerSubject = 100;
    
    %Call Function
    GetFullConnmat(inputPath, outputPath, subjRange, splitPerSubject);
end

%% Get FS Connectivity and Waytotal Matrices from full connectivity matrix above
if run == 2
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
if run == 3
    %Parameters
    inputMaskPath = '/net/store/nbp/projects/phasesim/databases/Bastian_DTI/masks/';
    inputcomplFSPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_overlappingFS/'];
    subjTotal = 22;
    
    %Call Function
    removeOverlappingFS(inputMaskPath, inputcomplFSPath, outputPath, subjTotal);
    
end
%% Get Recursive Ncut from full connectivity matrix
if run == 4
    
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

if run == 5
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

if run==6
    
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

if run == 7
    
    %CallFunction
    clusterStats();
end

%% Collect subject stats and metrics to one file

if run == 8
    filepath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150920StatsAndMetrics/';
    subjRange = 1:11;
    
    collectstatsandmetrics(filepath, subjRange);
end

%% Plot subjects stats and metrics

if run == 9
    
    filepath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150920StatsAndMetrics/';
    subjRange = 1:11;
    
    plotstatsandmetrics(filepath, subjRange);
end

%% Make clusternii

if run == 10
    
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

if run == 11
    for WeighingFactor = 0:0.1:1
        subjRange = [1];%:4, 6:13, 15, 17:22];%1
        decayParam = -1; % -0.25, -1
        %     WeighingFactor = 0; % 0.7
        normBy = 'sum';
        calcDistMat = 0;
        onlyChangeNorm = 1;
        InputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160417_Allconnmat';
        OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/decay'...
            num2str(decayParam) 'weigh' num2str(WeighingFactor)]  ;
        MainPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub'];
        useCosineSim = 0;
        
        CompSimMatCalcFn(subjRange, InputPath, OutputPath, MainPath, decayParam, WeighingFactor, calcDistMat, normBy, onlyChangeNorm, useCosineSim);  %% copied from DistanceMatrixCalc workdir/Arushi/20160217DistanceWeighingolddata
    end
end

%% Graclus recursive Ncut

if run == 12
    
    for WeighingFactor = 0:0.1:0.9
        for clusterCount = 1000
            disp('do clustering')
            recursiveSplit = true;
            %         clusterCount = 1000;
            subjRange = [1];%, 6:13, 15, 17:22];%1
            decayParam = -1; % -0.5, -1
            %     WeighingFactor = 0; % 0.7    if recursiveSplit
            normBy = 'sum';
            if recursiveSplit
                splitType = '';
            else
                splitType = 'NonRec';
            end
            InputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor)]; %'/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSubwithoutnormbymean/decay-0.25weigh0.5'; %
            OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor) 'conn' splitType];
            thresholdingFactor = 10;
            disp( ['WEIGHINGFACTOR:' num2str(WeighingFactor)]);
            graclusNcut(recursiveSplit, clusterCount, subjRange, InputPath, OutputPath, normBy, thresholdingFactor)
            
        end
    end
    
end

%% Global cut value calculation
if run == 13
    clustRange = 2:1000;
    subjRange = 1;
    compSimMatPathInit = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub';
    clusterIdPathInit = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut';
    decayParam = -1;
    WeighFactor = 0.5;
    recursiveSplit = true;
    useCosineSim = false;
    normBy = 'sum';
    thresholdFactor = 10;
    
    if useCosineSim
        text = 'cos';
    else
        text = 'conn';
    end
    
    if recursiveSplit
        recText = '';
    else
        recText = 'NonRec';
    end
    
    for subjNum = subjRange
        
        compSimMatPath = [compSimMatPathInit '/decay' num2str(decayParam) 'weigh' num2str(WeighFactor)...
            '/compSimMatWholeMax/compSimMatWholeMax' normBy num2str(subjNum)];
        GraclusPath = [clusterIdPathInit '/decay' num2str(decayParam) 'weigh' num2str(WeighFactor) text ...
            recText 'normby' normBy 'thresh' num2str(thresholdFactor)];
        clusterIdPath = [GraclusPath '/graclusResultnormBy' normBy 'thresh' num2str(thresholdFactor) 'subj' num2str(subjNum)];
        GlobalCutValCalc( clustRange, compSimMatPath, GraclusPath, clusterIdPath, thresholdFactor )
    end
end

%% Generate clusterwise connmat for all number of clusters
if run ==14
    subjRange = 1;
    recursiveSplit = true;
    decayParam = -1;
    WeighingFactor = 0.5;
    useCosineSim = false;
    normBy = 'sum';
    ThresholdFactor = 10;
    clusterRange = 2:1000;
    
    if useCosineSim
        conntext = 'cos';
    else
        conntext = 'conn';
    end
    if recursiveSplit
        recText = '';
    else
        recText = 'NonRec';
    end
    for subjNum = subjRange
        SPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/SSymDiag0/SSymDiag0subj' num2str(subjNum)];
        if recursiveSplit
        clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/decay'...
            num2str(decayParam) 'weigh' num2str(WeighingFactor) conntext recText 'normby' normBy 'thresh' num2str(ThresholdFactor)...
            '/graclusResultnormBy' normBy 'thresh' num2str(ThresholdFactor) 'subj' num2str(subjNum)];
        else
            clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160419_GraclusCut/decay'...
            num2str(decayParam) 'weigh' num2str(WeighingFactor) conntext recText 'normby' normBy 'thresh' num2str(ThresholdFactor)...
            '/gracluscollectedresult'];
        end
        CoordPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/FinalCoord/FinalCoord' num2str(subjNum)];
        OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/' recText 'Subj'  num2str(subjNum)];
        ZeroVoxelIdPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/zeroVoxelIdx/zeroBothVoxelIdx' num2str(subjNum)];
        GenerateClusterConnmat(clusterRange, clusterPath, CoordPath, SPath, recursiveSplit, OutputPath);
    end
end
%% Calculate cluster stats
if run==15
    subjRange = 1;
    clustRange = 2:1000;
    for subjNum = subjRange
        
        clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/details.mat'];
        
        clusterStatsNew(clusterPath, clustRange);
    end
end

%% Calculate intersecting clusters

if run == 16
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

if run == 17
    subjRange = 1;
    clustRange = 2:1000;
    for subjNum = subjRange
        for weighingFactor = 0:0.1:0.9
            
        for clusterNum = clustRange
            
            detClust = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/ClustOrgW'...
                num2str(WeighingFactor) '/detailsClust' num2str(clusterNum) '.mat']);
            number_of_outliers = zeros(clusterNum,1);
            number_of_connectedComponents = zeros(clusterNum,1);
            for singleCluster = 1:clusterNum
                voxel_coords = voxelCoordByCluster{singleCluster};
    [number_of_outliers(singleCluster), number_of_connectedComponents(singleCluster),clusterSecondMoment] = connectedComp(voxel_coords);
            end
        end
        end
    end
end
    
    %% Reduce stats and metrics to one file
   if run == 18
       clusterRange = 2:1000;
       
       for clusterNum = clusterRange
            StatTemp= load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/StatandMet/clusterstatandmetClust' num2str(clusterNum) '.mat']);
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
            save('/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj1/Rec/StatandMet/statsAll', 'statsAll', '-v7.3');
   end
   
   %% Plot Metrics
   if run == 19
   
   
  
  
  
   
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