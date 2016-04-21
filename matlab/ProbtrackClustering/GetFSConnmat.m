function  GetFSConnmat(inputPath,inputMaskPath, inputcomplFSPath, outputPath,subjTotal)

if ~inputPath
    inputPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150821_FullConnmat';
end

if ~inputMaskPath
    inputMaskPath = '/net/store/nbp/projects/phasesim/databases/Bastian_DTI/masks/';
end

if ~inputcomplFSPath
    inputcomplFSPath = '/net/store/nbp/projects/phasesim/workdir/Arushi/20150727Allsubjecttracking/';
end

if ~outputPath
    outputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/' num2str(yyyymmdd(datetime)) '_FSConnmat/'];
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
        connmat = load([inputPath '/connmatSubj' num2str(subjNum) '.mat']);
        connmat = connmat.connmat;
        fs_masks = cell(length(masks)-2,1);
        fs_masks_coord = cell(length(masks)-2,1);
        connmatfs = zeros(length(masks)-2,length(masks)-2);
       
        for n = 1:(length(masks)-2)
            fs_masks{n} = load_untouch_nii([finalMaskPath masks(n+2).name]);
            fs_masks_coord{n} = find(fs_masks{n}.img(find(complFSMask.img)));
        end
        for n = 1:length(masks)-2
            disp(n)
            k = connmat(fs_masks_coord{n},:);
    
            for m = 1:66
                j = k(:,fs_masks_coord{m});
                connmatfs(n,m) = sum(j(:));
            end
        end
        save([outputPath 'connmatFSsubj' caNum], 'connmatfs');
    end
end

end