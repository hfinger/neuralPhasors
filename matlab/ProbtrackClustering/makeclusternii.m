function makeclusternii(clusterfilepath, maskfilepath, savepath, modeRange, subjRange, clustRange)

for subjNum = subjRange
    for mode = modeRange
        switch mode
            case 1
                modeText = 'fsconn';
            case 2
                modeText = 'fscos';
            case 3
                modeText = 'fullconn';
            case 4
                modeText = 'fullcos';
        end
        if subjNum < 10
            caNum = ['0' num2str(subjNum)];
        else
            caNum = num2str(subjNum);
        end
        finalmaskpath = [maskfilepath 'ca' caNum 'FA_masks_FA_thr_012compl_fs_mask.nii'];
        usevoxelpath = [clusterfilepath num2str(subjNum) '/'];
        for clustNum = clustRange
            finalclusterpath = [clusterfilepath num2str(subjNum) '/' modeText '/' modeText 'Cluster' num2str(clustNum) '.mat'];
            mask = load_untouch_nii(finalmaskpath);
            clustermat = load(finalclusterpath);
            clustermat = clustermat.Cluster;
            
            for cluster = 1:clustNum
                tempVoxel = find(clustermat == cluster);
                n = 1;
                while 1
                    if exist([usevoxelpath '/useVoxelIdx' num2str(n) '.mat'], 'file')
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
                    
                    useVoxelIdx = load([usevoxelpath '/useVoxelIdx' num2str(n) '.mat']);
                    useVoxelIdx = useVoxelIdx.useVoxelIdx;
                    tempVoxel = useVoxelIdx(tempVoxel);
                    
                end
                
                cortex = find(mask.img);
                mask.img(cortex(tempVoxel)) = cluster*10;
                
               
                
            end
            mask.img((mask.img ==1)) = 0;
             mkdir([savepath num2str(subjNum) '/']);
             save_untouch_nii(mask, [savepath num2str(subjNum) '/' num2str(clustNum)]);
        end
    end
    
end
end
