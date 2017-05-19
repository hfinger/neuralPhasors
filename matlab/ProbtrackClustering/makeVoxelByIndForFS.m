inputMaskPath = '/net/store/nbp/projects/phasesim/databases/Bastian_DTI/masks/';
inputcomplFSPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
removeOverlapping = false;

outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_voxelByIndFS/'];

subjTotal = 22;

if ~exist(outputPath, 'dir')
    mkdir(outputPath);
end

for subjNum = 1:subjTotal
    if subjNum <10
        caNum = ['0' num2str(subjNum)];
    else
        caNum = num2str(subjNum);
    end
    
    finalMaskPath = [inputMaskPath 'ca' caNum '_1/FA_masks_FA_thr_012/'];
    if exist(finalMaskPath, 'dir')
        masks = dir(finalMaskPath);
        complFSMaskpath = [inputcomplFSPath 'ca' caNum 'FA_masks_FA_thr_012compl_fs_mask.nii'];
        complFSMask = load_untouch_nii(complFSMaskpath );
        zerovoxelIdx = '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/zeroVoxelIdx
        fsroi = zeros(sum(complFSMask.img(:))-length(zerovoxelIdx),1);
        voxelByInd = cell(length(masks)-2,1);
        for currentROI = 1:(length(masks)-2)
            currentROIMasks = load_untouch_nii([finalMaskPath masks(currentROI+2).name]);
            currentROIInd = find(currentROIMasks.img(find(complFSMask.img)));
            currentROISize = length(currentROIInd);
            if removeOverlapping
                for ROIInd = 1:length(currentROIInd)
                    if fsroi(currentROIInd(ROIInd))
                        size_other_fs = length(find(fsroi == fsroi(currentROIInd(ROIInd))));
                        if size_other_fs > currentROISize
                            fsroi(currentROIInd(ROIInd),1) = currentROI;
                        else
                            currentROISize = currentROISize -1;
                        end
                    else
                        fsroi(currentROIInd(ROIInd),1) = currentROI;
                    end
                end
                voxelByInd{currentROI} = find(fsroi == currentROI);
            else
                voxelByInd{currentROI} = currentROIInd;
            end
        end
    end
    if removeOverlapping
        save([outputPath caNum 'fsroiRemoveOverlap'],'fsroi');
        save([outputPath caNum 'voxelByIndFSRemoveOverlap'], 'voxelByInd');
    else
        save([outputPath caNum 'fsroiWithOverlap'],'fsroi');
        save([outputPath caNum 'voxelByIndFSWithOverlap'], 'voxelByInd');
    end
end