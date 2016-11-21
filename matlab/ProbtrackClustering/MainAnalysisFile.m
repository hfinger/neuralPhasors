% Generate connectivity matrices from tracking data for all subjects, cluster them and analyse them
% 20160407 Added CompSimMatCalc

%% Parameters
run = 12;
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
    
    for WeighingFactor = 0:0.1:1
        if WeighingFactor == 0.5
            continue;
        end
        disp('do clustering')
        recursiveSplit = false;
        clusterCount = 1000;
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