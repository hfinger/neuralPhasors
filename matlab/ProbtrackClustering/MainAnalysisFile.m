%% Generate connectivity matrices from tracking data for all subjects, cluster them and analyse them

%% Get Full Connectivity and Waytotal matrices from grid job results of ProbtrackX
%Parameters
inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150804_ProbtrackXallsubjects';
outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_FullConnmat'];
subjTotal = 22;
splitPerSubject = 50;

%Call Function
GetFullConnmat(inputPath, outputPath, subjTotal, splitPerSubject);

%% Get FS Connectivity and Waytotal Matrices from full connectivity matrix above


