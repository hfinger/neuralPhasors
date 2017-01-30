function GenerateClusterConnmat(clusterRange, clusterPath, CoordPath, SPath, recursiveSplit, OutputPath)
%Modified from RecursiveNcutPostProcess
%   Detailed explanation goes here

S = load(SPath);
S = S.S;

% zeroBothVoxelIdx = load(ZeroVoxelIdPath);
% zeroBothVoxelIdx = zeroBothVoxelIdx.zeroBothVoxelIdx;
% 
% S(zeroBothVoxelIdx,:) = [];
% S(:,zeroBothVoxelIdx) = [];

FinalCoord = load(CoordPath);
FinalCoord = FinalCoord.FinalCoord;
Cluster = load(clusterPath);
stepWiseClusters = Cluster.stepWiseClusters;

clusterCenter = cell(clusterRange(end));
clusterConnmat = cell(clusterRange(end),1);
normClusterConnmat = cell(clusterRange(end),1);
voxelIndByCluster = cell(clusterRange(end));
voxelCoordByCluster = cell(clusterRange(end));

for clusterNum = clusterRange
    clustConnmat = zeros(clusterNum);
    normClustConnmat = zeros(clusterNum);
    
    if recursiveSplit
    if clusterNum~=2
        clustConnmat(1:clusterNum-1,1:clusterNum-1) = clusterConnmat{clusterNum-1};
        normClustConnmat(1:clusterNum-1,1:clusterNum-1) = normClusterConnmat{clusterNum-1};
        voxelIndByCluster(:,clusterNum) = voxelIndByCluster(:,clusterNum-1);
        clusterCenter(:,clusterNum) = clusterCenter(:, clusterNum-1);
        voxelCoordByCluster(:,clusterNum) = voxelCoordByCluster(:,clusterNum-1);
    end
    newClust = stepWiseClusters(:,clusterNum) - stepWiseClusters(:,clusterNum-1);
    newClustInd = find(newClust ~=0);
    
    OldClustNum = stepWiseClusters(newClustInd(1),(clusterNum-1));
    NewClustNum = stepWiseClusters(newClustInd(1), clusterNum);
    
    voxelsInOldClust = find(stepWiseClusters(:,clusterNum) == OldClustNum);
    voxelCoordByCluster{OldClustNum,clusterNum}  = FinalCoord(voxelsInOldClust,:);
    clusterCenter{OldClustNum,clusterNum} = mean(voxelCoordByCluster{OldClustNum,clusterNum});
    voxelIndByCluster{OldClustNum,clusterNum} = voxelsInOldClust;
    
    voxelsInNewClust = find(stepWiseClusters(:,clusterNum) == NewClustNum);
    voxelCoordByCluster{NewClustNum,clusterNum}  = FinalCoord(voxelsInNewClust,:);
    clusterCenter{NewClustNum,clusterNum} = mean(voxelCoordByCluster{NewClustNum,clusterNum});
    voxelIndByCluster{NewClustNum,clusterNum} = voxelsInNewClust;
    
    interConnmatOld1 = sum(S(voxelsInOldClust,:),1);
    interConnmatNew1 = sum(S(voxelsInNewClust,:),1);
    
    interConnmatOld2 = sum(S(:,voxelsInOldClust),2);
    interConnmatNew2 = sum(S(:,voxelsInNewClust),2);
        
    for singleCluster = 1:clusterNum
        clustConnmat(singleCluster,OldClustNum) = sum(interConnmatOld2(voxelIndByCluster{singleCluster,clusterNum}));
        clustConnmat(singleCluster,NewClustNum) = sum(interConnmatNew2(voxelIndByCluster{singleCluster,clusterNum}));
        clustConnmat(OldClustNum,singleCluster) = sum(interConnmatOld1(voxelIndByCluster{singleCluster,clusterNum}));
        clustConnmat(NewClustNum,singleCluster) = sum(interConnmatNew1(voxelIndByCluster{singleCluster,clusterNum}));
        
        normClustConnmat(singleCluster,OldClustNum) = sum(interConnmatOld2(voxelIndByCluster{singleCluster,clusterNum}))/...
            (size(voxelIndByCluster{singleCluster,clusterNum},1)*(size(voxelIndByCluster{OldClustNum,clusterNum},1)));
        normClustConnmat(singleCluster,NewClustNum) = sum(interConnmatNew2(voxelIndByCluster{singleCluster,clusterNum}))/...
            (size(voxelIndByCluster{singleCluster,clusterNum},1)*(size(voxelIndByCluster{NewClustNum,clusterNum},1)));
         normClustConnmat(OldClustNum,singleCluster) = sum(interConnmatOld2(voxelIndByCluster{singleCluster,clusterNum}))/...
            (size(voxelIndByCluster{singleCluster,clusterNum},1)*(size(voxelIndByCluster{OldClustNum,clusterNum},1)));
        normClustConnmat(NewClustNum,singleCluster) = sum(interConnmatNew2(voxelIndByCluster{singleCluster,clusterNum}))/...
            (size(voxelIndByCluster{singleCluster,clusterNum},1)*(size(voxelIndByCluster{NewClustNum,clusterNum},1)));
        
    end
       clusterConnmat{clusterNum} = clustConnmat;
        normClusterConnmat{clusterNum} = normClustConnmat;   
    else
            voxelIndByCluster{clusterNum} = arrayfun(@(x) find(stepWiseClusters{clusterNum} == x), 1:clusterNum, 'UniformOutput', false)';
            voxelCoordByCluster{clusterNum} = cellfun(@(x) FinalCoord(x,:), voxelIndByCluster{clusterNum}, 'UniformOutput', false);
            clusterCenter{clusterNum} = cellfun(@(x) mean(x,1), voxelCoordByCluster{clusterNum}, 'UniformOutput', false);
            disp(num2str(clusterNum));
        for singleCluster = 1:clusterNum
            
            for otherCluster = 1:clusterNum
                a = S(voxelIndByCluster{clusterNum}{singleCluster},voxelIndByCluster{clusterNum}{otherCluster});
                clustConnmat(singleCluster, otherCluster) = sum(a(:));
                normClustConnmat(singleCluster, otherCluster) = clustConnmat(singleCluster, otherCluster)/(size(voxelIndByCluster{clusterNum}{singleCluster},1) *size(voxelIndByCluster{clusterNum}{otherCluster},1));
            end
        end
        
        clusterConnmat{clusterNum} = clustConnmat;
        normClusterConnmat{clusterNum} = normClustConnmat;   

            
            
    end
  
end

if ~exist(OutputPath, 'dir');
    mkdir(OutputPath);
end
save([OutputPath '/details' ],  'clusterCenter','clusterConnmat',...
    'normClusterConnmat', 'voxelIndByCluster', 'voxelCoordByCluster', '-v7.3')

end

