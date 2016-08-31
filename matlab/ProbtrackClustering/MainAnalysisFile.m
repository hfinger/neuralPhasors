% Generate connectivity matrices from tracking data for all subjects, cluster them and analyse them
% 20160407 Added CompSimMatCalc

%% Parameters
run = 11;
%% Get Full Connectivity and Waytotal matrices from grid job results of ProbtrackX
if run == 1
    %Parameters
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150804_ProbtrackXallsubjects';
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_FullConnmat'];
    subjTotal = 22;
    splitPerSubject = 50;
    
    %Call Function
    GetFullConnmat(inputPath, outputPath, subjTotal, splitPerSubject);
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
    subjRange = [1:4, 6:13, 15, 17:22];%1
    decayParam = -0.5; % -0.25, -1
    WeighingFactor = 0.7; % 0.7
    normBy = 'sum';
    calcDistMat = 1;
    onlyChangeNorm = 0;
    InputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSub';
    OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSub/decay'...
        num2str(decayParam) 'weigh' num2str(WeighingFactor)]  ;
    
    CompSimMatCalc(subjRange, InputPath, OutputPath, decayParam, WeighingFactor, calcDistMat, normBy, onlyChangeNorm);  %% copied from DistanceMatrixCalc workdir/Arushi/20160217DistanceWeighingolddata
    
end

%% Graclus recursive Ncut

if run == 12
    disp('do clustering')
    recursiveSplit = true;
    clusterCount = 1000;
    subjRange = [1:4, 6:13, 15, 17:22];%, 6:13, 15, 17:22];%1
    decayParam = -0.25; % -0.5, -1
    WeighingFactor = 0.5; % 0.7    if recursiveSplit
    normBy = 'sum';
    InputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSub/decay'...
        num2str(decayParam) 'weigh' num2str(WeighingFactor)]; %'/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSubwithoutnormbymean/decay-0.25weigh0.5'; %
    OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160411_GraclusCut/decay'...
        num2str(decayParam) 'weigh' num2str(WeighingFactor) 'conn'];
    
    graclusNcut(recursiveSplit, clusterCount, subjRange, InputPath, OutputPath, normBy)
    
    
    
end