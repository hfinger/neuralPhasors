function RecursiveNcutPostProcess(subjTotal, clusterTotal, clusterPath, clusterType, maskPath, connmatPath)
%RECURSIVENCUTPOSTPROCESS Summary of this function goes here
%   Detailed explanation goes here

if ~clusterPath
    clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150824Allsubjectrecursivencut/'];
end

if ~clusterTotal
    clusterTotal = 1000;
end

switch clusterType
    case 'fsconn'
        clusterStart = 67;
    case 'fscos'
        clusterStart = 67;
    case 'fullconn'
        clusterStart = 2;
    case 'fullcos'
        clusterStart = 2;
end

if ~maskPath
    maskPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
end

if~connmatPath
    connmatPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20150821_FullConnmat'];
end


for subjNum = 1:subjTotal
    folder = [clusterPath num2str(subjNum) '/' clusterType '/'];
    if subjNum <10
        caNum = ['0' num2str(subjNum)];
    else
        caNum = num2str(subjNum);
    end
    complFSMaskpath = [maskPath 'ca' caNum 'FA_masks_FA_thr_012compl_fs_mask.nii'];
    complFSMask = load_untouch_nii(complFSMaskpath );
    [complCoordx, complCoordy, complCoordz] = ind2sub(size(complFSMask.img),find(complFSMask.img));
    complCoords = [complCoordx complCoordy complCoordz];
    if ~exist([connmatPath '/connmatSubj' num2str(subjNum) '.mat'], 'file')
        continue
    end
    connmat = load([connmatPath '/connmatSubj' num2str(subjNum) '.mat']);
    connmat = connmat.connmat;
    
    for clusterNum = clusterStart:clusterTotal+1
        disp(['Subj' num2str(subjNum) 'cluster' num2str(clusterNum)]);
        Cluster = load([folder clusterType 'Cluster' num2str(clusterNum) '.mat']);
        Cluster = Cluster.Cluster;
        clusterCenter = zeros(clusterNum,3);
        clusterConnmat = zeros(clusterNum, clusterNum);
        voxeltarget = cell(clusterNum,1);
        flow = 1;
        n = 1;
        while 1
            if exist([clusterPath num2str(subjNum) '/useVoxelIdx' num2str(n) '.mat'], 'file')
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
            
            useVoxelIdx = load([clusterPath num2str(subjNum) '/useVoxelIdx' num2str(n) '.mat']);
            useVoxelIdx = useVoxelIdx.useVoxelIdx;
            if flow
                tempVoxel = useVoxelIdx;
                flow = 0;
            else
                tempVoxel = useVoxelIdx(tempVoxel);
            end
        end
        
        for singleCluster = 1:clusterNum
            voxeltarget{singleCluster} = tempVoxel(find(Cluster==singleCluster));
        end
        
        for singleCluster = 1:clusterNum
            currentCluster = complCoords(voxeltarget{singleCluster},1:3);
            clusterCenter(singleCluster,1:3) = mean(currentCluster,1);
            tempconnmat = sum(connmat(voxeltarget{singleCluster},:),1);
            for otherCluster = 1:clusterNum
                clusterConnmat(singleCluster,otherCluster) = sum(tempconnmat(:,voxeltarget{otherCluster}));
            end
            %      if ~(mod(clusterNum, 100))
            %             compl_fs_mask.img(currentCluster) = singleCluster*30;
            %             save_untouch_nii(compl_fs_mask,'compl_fs_mask');
            %      end
        end
        if ~exist([folder clusterType '/postprocessing'], 'dir');
            mkdir([folder clusterType '/postprocessing']);
        end
        save([folder clusterType '/postprocessing/clusterConnmat',num2str(clusterNum)], 'clusterConnmat');
        save([folder clusterType '/postprocessing/clusterCenter' num2str(clusterNum)],'clusterCenter');
        clear clusterCenter;
    end
end
end

