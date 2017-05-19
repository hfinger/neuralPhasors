subjRange = [1];
decayParam = -1;
useCosineSim = false;
normBy = 'sum';
WholeNormText = 'WholeMax';
WeighFacRange = [0, 0.5, 1];
RecursiveText = 'Rec';
threshRange = [100,10,5,1];
clustRange = [2:8];
FinalCoordPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/FinalCoord/';
outputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160706_ClusterNii/';


if useCosineSim
    cosText = 'cos';
else
    cosText = 'conn';
end
for subjNum = subjRange
    caNum = ['0' num2str(subjNum)];
    complFSMaskfilename = ['ca' num2str(caNum) 'FA_masks_FA_thr_012compl_fs_mask'];
    complFSMask = load_untouch_nii(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160407_CompSimMatCalcNewSub/complFSMask/' complFSMaskfilename '.nii']);
    FinalCoord = load([FinalCoordPath 'FinalCoord' num2str(subjNum) '.mat']);
    FinalCoord = FinalCoord.FinalCoord;
    for WeighFacNum = WeighFacRange
%         clusterfilepath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/'...
%             RecursiveText '/decay' num2str(decayParam) 'weigh' num2str(WeighFacNum) ...
%             '/normby' normBy 'thresh'];
        clusterfilepath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160706_GraclusCuttest/conn/Rec/decay-1weigh0/1/normbysumthresh'];
        for threshFactor = threshRange
            if exist([clusterfilepath num2str(threshFactor) ], 'dir')
%                            if exist([clusterfilepath num2str(threshFactor) '/' cosText '/Subj' num2str(subjNum)], 'dir')

                threshFound = 1;
                break;
            end
        end
        
        if ~threshFound
            error('No cluster matrix found for subj %i for weighingFactor %i',subjNum, WeighingFactor);
        end
                clusterfilepath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160706_ClusteringPostprocessingtest/Rec/decay-1weigh0/normbysumthresh10/conn/Subj1' '/detailsClust2to8'  '.mat'];

%         clusterfilepath = [clusterfilepath num2str(threshFactor) '/' cosText '/Subj' num2str(subjNum) '/detailsClust'  '.mat'];
            % num2str(clustRange(1)) 'to' num2str(clustRange(end))
            
        cluster = load(clusterfilepath);
        
        for clustNum = clustRange
            clusterCoords = cluster.voxelCoordByCluster(:,clustNum);
            complFSMask.img(:) = 0;
            for currentClusterNum = 1:clustNum
                currentCoords = clusterCoords{currentClusterNum};
                currentInd = sub2ind(size(complFSMask.img),currentCoords(:,1), currentCoords(:,2), currentCoords(:,3));
                complFSMask.img(currentInd) = currentClusterNum;
            end
            FinalOutputPath = [outputPath RecursiveText '/decay' num2str(decayParam) 'weigh' num2str(WeighFacNum) ...
                '/normby' normBy 'thresh' num2str(threshFactor) '/' ...
                cosText '/Subj' num2str(subjNum) '/'];
            if ~exist(FinalOutputPath, 'dir');
                mkdir(FinalOutputPath);
            end
            save_untouch_nii(complFSMask, [FinalOutputPath 'cluster' num2str(clustNum)]);
        end
    end
    
end
