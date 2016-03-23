function GetRecursiveNcut(inputPath, clusterPath, useCosineSimilarity, useFSROI, numberIterations, outputPath, subjTotal )
%GetRecursiveNcut Perform Ncut for stated iterations for all subjects
%   modified version of RecursiveNcut file to make it suitable for all
%   subjects

%% parameters
if ~numberIterations 
    numberIterations = 1000;
end
if ~inputPath
        inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150821_FullConnmat';
end

if  ~outputPath  
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/'];
end

if ~clusterPath
    clusterPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150824_overlappingFS/';
end

%% Perform function
for subjNum = 19:subjTotal
    
    if subjNum <10
        caNum = ['0' num2str(subjNum)];
    else
        caNum = num2str(subjNum);
    end
% load matrix
connmatpath = [inputPath '/connmatSubj' num2str(subjNum) '.mat'];
if exist(connmatpath, 'file');
    fullconnmat =  load(connmatpath);
else
    continue;
end
fullconnmat = fullconnmat.connmat;

%% Remove zeros
previousSize = 0;
n = 1;
while 1
    useVoxelIdx = find(sum(fullconnmat,2)~=0);
    fullconnmat = fullconnmat(useVoxelIdx,useVoxelIdx);
    if size(fullconnmat,1) == previousSize
        break;
    else
        previousSize = size(fullconnmat,1);
        save ([outputPath num2str(subjNum) '/useVoxelIdx' num2str(n) '.mat'], 'useVoxelIdx');
        n = n + 1;
    end
end
    

%% Normalise by sum of rows
fullconnmat = bsxfun(@rdivide,fullconnmat,sum(fullconnmat,2));

%% Make symmetric
fullconnmat = fullconnmat + fullconnmat';

%% Apply cosine similarity if it is the case
if useCosineSimilarity
    tmp = fullconnmat*fullconnmat';
    tmp = bsxfun(@rdivide,tmp, sqrt(sum(fullconnmat.^2,2)) );
    fullconnmat = bsxfun(@rdivide,tmp, sqrt(sum(fullconnmat.^2,1)) );
    clear tmp;
end
%% Make diagonal 1
setDiagTo = 1;
fullconnmat = fullconnmat-diag(diag(fullconnmat));
b = sparse(diag(diag(fullconnmat) + setDiagTo));
fullconnmat = fullconnmat + b;
clear b;


%% Final matrix

    final  = full(fullconnmat);


%% whole matrix - connectivity
if useFSROI == 0 && useCosineSimilarity == 0
    startClusterNum = 1;
    Cluster = ones(size(final,1),1);
    clear fullconnmat;
    text = 'fullconn';
end

%% whole matrix - cosine similarity
if useFSROI == 0 && useCosineSimilarity
    startClusterNum = 1;
    Cluster = ones(size(final,1),1);
    clear cosS;
    clear fullconnmat;
    text = 'fullcos';
end


%% FS ROI
if useFSROI
    Cluster = load([clusterPath caNum 'fsroi.mat']);
    Cluster = Cluster.fsroi;
    n = 1;
    flow = 1;
    while 1
        if exist([outputPath num2str(subjNum) '/useVoxelIdx' num2str(n) '.mat'], 'file')
            n = n+1;
            continue;
        end
        break;
    end
    while 1
        n = n-1;
        if ~n
            break;
        end
        
        useVoxelIdx = load([outputPath num2str(subjNum) '/useVoxelIdx' num2str(n) '.mat']);
        useVoxelIdx = useVoxelIdx.useVoxelIdx;
        if flow
            tempVoxel = useVoxelIdx;
            flow = 0;
        else
            tempVoxel = useVoxelIdx(tempVoxel);
        end
    end
    
    Cluster = Cluster(tempVoxel);
end

%% 66 ROI - connectivity
if useFSROI && useCosineSimilarity == 0
    startClusterNum = 66;
    clear fullconnmat;
    text = 'fsconn';
end

%% 66 ROI - cosine similarity
if useFSROI && useCosineSimilarity
    startClusterNum = 66;
    clear cosS;
    clear fullconnmat;
    text = 'fscos';
end

%% Recursive ncut

for iteration = startClusterNum:numberIterations
    greatestSum = 0;
    greatestRoi = 0;
    for j = 1:iteration
        currentSum = size(find(Cluster == j),1);
        if currentSum > greatestSum
            greatestSum = currentSum;
            greatestRoi = j;
        end
    end
    currentRoiVoxel = find(Cluster == greatestRoi);
    if size(Cluster,1) == size(currentRoiVoxel,1)
        currentFinal = final;
    else
        currentFinal = final(currentRoiVoxel, currentRoiVoxel);
      
    end
      
    nbCluster = 2;
    [NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(currentFinal,nbCluster);
    Cluster(currentRoiVoxel(NcutDiscrete(:,1) == 0)) = iteration+1;
    clear NcutDiscrete;
    clear NcutEigenvectors;
    clear NcutEigenvalues;
    clear currentFinal;
    folder = [outputPath num2str(subjNum) '/' text '/'];
    if ~exist(folder, 'dir');
        mkdir(folder);
    end
    save([ folder text 'Cluster' num2str(iteration+1)], 'Cluster');

end

clear fullconnmat;
    



end
end

