%% Generate connectivity matrices from tracking data for all subjects, cluster them and analyse them
run = 5;
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