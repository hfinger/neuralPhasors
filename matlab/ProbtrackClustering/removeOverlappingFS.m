function removeOverlappingFS( inputMaskPath, inputcomplFSPath, outputPath,subjTotal )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ~inputMaskPath
    inputMaskPath = '/net/store/nbp/projects/phasesim/databases/Bastian_DTI/masks/';
end

if ~inputcomplFSPath
    inputcomplFSPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
end

if ~outputPath
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_overlappingFS/'];
end

if ~subjTotal
    subjTotal = 22;
end

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
        fsroi = zeros(sum(complFSMask.img(:)),1);
        for currentROI = 1:(length(masks)-2)
            currentROIMasks = load_untouch_nii([finalMaskPath masks(currentROI+2).name]);
            currentROICoord = find(currentROIMasks.img(find(complFSMask.img)));
            currentROISize = length(currentROICoord);
            for ROICoord = 1:length(currentROICoord)
                if fsroi(currentROICoord(ROICoord))
                    size_other_fs = length(find(fsroi == fsroi(currentROICoord(ROICoord))));
                    if size_other_fs > currentROISize
                        fsroi(currentROICoord(ROICoord),1) = currentROI;
                    else
                        currentROISize = currentROISize -1;
                    end
                else
                    fsroi(currentROICoord(ROICoord),1) = currentROI;
                end
            end
        end
        save([outputPath caNum 'fsroi'],'fsroi');
    end
end
end
