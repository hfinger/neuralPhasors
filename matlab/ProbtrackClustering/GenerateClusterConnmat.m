function GenerateClusterConnmat(splitRange, clusterPath, CoordPath, SPath, recursiveSplit, OutputPath)
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

clusterCenter = cell(splitRange(end));
clusterConnmat = cell(splitRange(end),1);
normClusterConnmat = cell(splitRange(end),1);
voxelIndByCluster = cell(splitRange(end));
voxelCoordByCluster = cell(splitRange(end));

for split = splitRange
    disp(num2str(split));
    clustConnmat = zeros(split);
    normClustConnmat = zeros(split);
    
    if recursiveSplit
        if split~=2
            clustConnmat(1:split-1,1:split-1) = clusterConnmat{split-1};
            normClustConnmat(1:split-1,1:split-1) = normClusterConnmat{split-1};
            voxelIndByCluster(:,split) = voxelIndByCluster(:,split-1);
            clusterCenter(:,split) = clusterCenter(:, split-1);
            voxelCoordByCluster(:,split) = voxelCoordByCluster(:,split-1);
        end
        newClust = stepWiseClusters(:,split) - stepWiseClusters(:,split-1);
        newClustInd = find(newClust);
        
        OldClustNum = stepWiseClusters(newClustInd(1),(split-1));
        NewClustNum = stepWiseClusters(newClustInd(1), split);
        
        voxelsInOldClust = find(stepWiseClusters(:,split) == OldClustNum);
        voxelCoordByCluster{OldClustNum,split}  = FinalCoord(voxelsInOldClust,:);
        clusterCenter{OldClustNum,split} = mean(voxelCoordByCluster{OldClustNum,split});
        voxelIndByCluster{OldClustNum,split} = voxelsInOldClust;
        
        voxelsInNewClust = find(stepWiseClusters(:,split) == NewClustNum);
        voxelCoordByCluster{NewClustNum,split}  = FinalCoord(voxelsInNewClust,:);
        clusterCenter{NewClustNum,split} = mean(voxelCoordByCluster{NewClustNum,split});
        voxelIndByCluster{NewClustNum,split} = voxelsInNewClust;
        
        interConnmatOld1 = sum(S(voxelsInOldClust,:),1);
        interConnmatNew1 = sum(S(voxelsInNewClust,:),1);
        
        interConnmatOld2 = sum(S(:,voxelsInOldClust),2);
        interConnmatNew2 = sum(S(:,voxelsInNewClust),2);
        
        for singleCluster = 1:split
            clustConnmat(singleCluster,OldClustNum) = sum(interConnmatOld2(voxelIndByCluster{singleCluster,split}));
            clustConnmat(singleCluster,NewClustNum) = sum(interConnmatNew2(voxelIndByCluster{singleCluster,split}));
            clustConnmat(OldClustNum,singleCluster) = sum(interConnmatOld1(voxelIndByCluster{singleCluster,split}));
            clustConnmat(NewClustNum,singleCluster) = sum(interConnmatNew1(voxelIndByCluster{singleCluster,split}));
            
            normClustConnmat(singleCluster,OldClustNum) = sum(interConnmatOld2(voxelIndByCluster{singleCluster,split}))/...
                (size(voxelIndByCluster{singleCluster,split},1)*(size(voxelIndByCluster{OldClustNum,split},1)));
            normClustConnmat(singleCluster,NewClustNum) = sum(interConnmatNew2(voxelIndByCluster{singleCluster,split}))/...
                (size(voxelIndByCluster{singleCluster,split},1)*(size(voxelIndByCluster{NewClustNum,split},1)));
            normClustConnmat(OldClustNum,singleCluster) = sum(interConnmatOld2(voxelIndByCluster{singleCluster,split}))/...
                (size(voxelIndByCluster{singleCluster,split},1)*(size(voxelIndByCluster{OldClustNum,split},1)));
            normClustConnmat(NewClustNum,singleCluster) = sum(interConnmatNew2(voxelIndByCluster{singleCluster,split}))/...
                (size(voxelIndByCluster{singleCluster,split},1)*(size(voxelIndByCluster{NewClustNum,split},1)));
            
        end
        clusterConnmat{split} = clustConnmat;
        normClusterConnmat{split} = normClustConnmat;
        
        
    else
        voxelIndByCluster(1:split,split) = arrayfun(@(x) find(stepWiseClusters(:,split) == x), 1:split, 'UniformOutput', false)';
        voxelCoordByCluster(1:split,split) = cellfun(@(x) FinalCoord(x,:), voxelIndByCluster(1:split,split), 'UniformOutput', false);
        clusterCenter(1:split,split) = cellfun(@(x) mean(x,1), voxelCoordByCluster(1:split,split), 'UniformOutput', false);
        disp(num2str(split));
        for singleCluster = 1:split
            
            for otherCluster = 1:split
                a = S(voxelIndByCluster{singleCluster, split},voxelIndByCluster{otherCluster,split});
                clustConnmat(singleCluster, otherCluster) = sum(a(:));
                normClustConnmat(singleCluster, otherCluster) = clustConnmat(singleCluster, otherCluster)/(size(voxelIndByCluster{singleCluster,split},1) *size(voxelIndByCluster{otherCluster,split},1));
            end
        end
        
        clusterConnmat{split} = clustConnmat;
        normClusterConnmat{split} = normClustConnmat;
        
        
   end
    
end

if ~exist(OutputPath , 'dir')
    mkdir(OutputPath);
end

save([OutputPath '/detailsClust' num2str(splitRange(1)) 'to' num2str(splitRange(end))] , 'clusterConnmat', 'normClusterConnmat', 'clusterCenter', 'voxelIndByCluster', 'voxelCoordByCluster', '-v7.3');
end





