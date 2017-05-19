prompt = {'Compare to FS?'...
    'Compare to next cluster?',...
    'Compare across different cosType?',...
    'Compare for different weighing Factors?',...
    'Compare for different clustering types(Rec or NonRec)?'};
dlg_title = 'Input only one true value';
num_lines = 1;
defaultans = {'true',...
    'false',...
    'false',...
    'false', ...
    'false'};


answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
    FScompare = str2num(answer{1});
    NextClusterCompare     = str2num(answer{2});
    CosTypeCompare    = str2num(answer{3});
    WeighingFactorCompare    = str2num(answer{4});
    ClustTypeCompare  = str2num(answer{5});
end

if FScompare
    prompt = {'Enter FS voxelByInd path:',...
        'Do you want to keep overlap(True or False)?',...
        'Enter Weighing Factor range:',...
        'Enter subject range:',...
        'Enter clustering range:',...
        'Enter Decay Parameter:',...
        'Do you want to use Cosine Similarity("true" or "false")?',...
        'Calculate for Recursive Split("true" or "false")?',...
        'Enter individual normalisation factor("sum" or "mean")',...
        'Enter whole normalisation("WholeMax")',...
        'Enter Thresholding Factor Search Range',...
        'Enter PostProcessing Path:',...
        'Output Path:'};
    dlg_title = 'FS Path';
    num_lines = 1;
    defaultans = {'/net/store/nbp/projects/phasesim/workdir/Arushi/20160606_voxelByIndFS/',...
        'false',...
        '(0:0.1:1)',...
        '[1]',...
        '(66:100)',...
        '-1', ...
        'false',...
        'true',...
        'sum',...
        'WholeMax',...
        '[100,10,5,1]',...
        '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_ClusteringPostprocessing/',...
        '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/201600709_IntersectingClusters/'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    FSPath        = answer{1};
    OverlapFlag   = str2num(answer{2});
    WeighFacRange = str2num(answer{3});
    subjRange     = str2num(answer{4});
    clustRange    = str2num(answer{5});
    decayParam    = answer{6};
    useCosineSim  = str2num(answer{7});
    recursiveSplit = str2num(answer{8});
    normBy         = answer{9};
    WholeNormText  = answer{10};
    threshRange    = str2num(answer{11});
    PostProcessPath = answer{12};
    OutputPath     = answer{13};
    
    if recursiveSplit
        RecText = 'Rec';
    else
        RecText = 'NonRec';
    end
    
    if useCosineSim
        cosText = 'cos';
    else
        cosText = 'conn';
    end
    
    for subjNum = subjRange
        if subjNum >= 10
            caNum = [num2str(subjNum)];
        else
            caNum = ['0' num2str(subjNum)];
        end
        
        if OverlapFlag
            OverlapText = 'WithOverlap';
        else
            OverlapText = 'RemoveOverlap';
        end
        FinalFSPath = [FSPath caNum 'voxelByIndFS' OverlapText '.mat'];
        voxelByInd = load(FinalFSPath);
        voxelByInd = voxelByInd.voxelByInd;
        
        
        for WeighFactor = WeighFacRange
            
            threshCount = 1;
            while 1
                
                FinalPostProcessPath = [PostProcessPath RecText '/decay' decayParam 'weigh'...
                    num2str(WeighFactor) '/normby' normBy 'thresh' ...
                    num2str(threshRange(threshCount)) '/' cosText '/Subj' num2str(subjNum) '/detailsClust.mat'];
                
                if ~exist(FinalPostProcessPath, 'file')
                    threshCount = threshCount + 1;
                    if threshCount > length(threshRange)
                        error(['No post processing file found for subj' num2str(subjNum)...
                            'with weighingFactor' num2str(WeighFactor) ...
                            'Recursive Type ' RecText...
                            'Cos Type ' cosText...
                            'normalised by' normBy]);
                    end
                else
                    break;
                end
            end
            
            detailsClust = load(FinalPostProcessPath);
            voxelIndAll = detailsClust.voxelIndByCluster;
            
            for clusterNum = clustRange
                
                voxelIndByCluster1 = voxelIndAll(:,clusterNum);
                voxelIndByCluster2 = voxelByInd;
                
                ret = calcCommonClusters(voxelIndByCluster1, voxelIndByCluster2);
                
                finalret.IntersectVal{clusterNum} =  ret.IntersectVal;
                finalret.IntersectInd{clusterNum,1} = ret.IntersectInd1;
                finalret.IntersectInd{clusterNum,2} =  ret.IntersectInd2;
                finalret.IntersectCount{clusterNum} = ret.IntersectCount;
                finalret.IntersectRatio{clusterNum,1} = ret.IntersectRatio1;
                finalret.IntersectRatio{clusterNum,2} = ret.IntersectRatio2;
                finalret.MaxIntRatio{clusterNum,1} = ret.MaxIntRatio1;
                finalret.MaxIntRatio{clusterNum,2} = ret.MaxIntRatio2;
                finalret.IndMaxIntRatio{clusterNum,1} = ret.IndMaxIntRatio1;
                finalret.IndMaxIntRatio{clusterNum,2} = ret.IndMaxIntRatio2;
                
                
            end
            ComparisonText = 'FSCompare';
            FinalOutputPath = [OutputPath ComparisonText OverlapText '/Subj' num2str(subjNum)... 
                '/' cosText '/' RecText '/decay' decayParam 'weigh' num2str(WeighFactor) '/'];
            
            if ~exist(FinalOutputPath, 'dir')
                mkdir(FinalOutputPath);
            end
            save([FinalOutputPath 'IntersectRange' num2str(clustRange(1)) 'to' num2str(clustRange(end))],... 
                  'finalret', '-v7.3');
              disp(['saved subj' num2str(subjNum) 'weigh' num2str(WeighFactor)]);
            
            
        end
    end
    
elseif NextClusterCompare
    prompt = {'Enter Weighing Factor range:',...
        'Enter subject range:',...
        'Enter clustering range:',...
        'Enter Decay Parameter:',...
        'Do you want to use Cosine Similarity("true" or "false")?',...
        'Calculate for Recursive Split("true" or "false")?',...
        'Enter individual normalisation factor("sum" or "mean")',...
        'Enter whole normalisation("WholeMax")',...
        'Enter Thresholding Factor Search Range'...
        'Enter PostProcessing Path:'};
    dlg_title = 'FS Path';
    num_lines = 1;
    defaultans = {'/net/store/nbp/projects/phasesim/workdir/Arushi/20160606_voxelByIndFS/',...
        'false',...
        '(0:0.1:1)',...
        '[1]',...
        '(66:100)',...
        '-1', ...
        'false',...
        'true',...
        'sum',...
        'WholeMax',...
        '[100,10,5,1]',...
        '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    WeighFacRange = str2num(answer{1});
    subjRange     = str2num(answer{2});
    clustRange    = str2num(answer{3});
    decayParam    = answer{4};
    useCosineSim  = str2num(answer{5});
    recursiveSplit = str2num(answer{6});
    normBy         = answer{7};
    WholeNormText  = answer{8};
    threshRange    = str2num(answer{9});
    PostProcessPath = answer{10};
    
    if recursiveSplit
        RecText = 'Rec';
    else
        RecText = 'NonRec';
    end
    
    if useCosineSim
        cosText = 'cos';
    else
        cosText = 'conn';
    end
    
    for subjNum = subjRange
        if subjNum >= 10
            caNum = [num2str(subjNum)];
        else
            caNum = ['0' num2str(subjNum)];
        end
        for WeighFactor = WeighFacRange
            
            threshCount = 1;
            while 1
                
                FinalPostProcessPath = [PostProcessPath RecText '/decay' decayParam 'weigh'...
                    num2str(WeighFactor) '/normby' normBy 'thresh' ...
                    num2str(threshRange(threshCount)) '/' cosText '/Subj' num2str(subjNum) '/detailsClust.mat'];
                
                if ~exist(FinalPostProcessPath, 'file')
                    threshCount = threshCount + 1;
                    if threshCount > length(threshRange)
                        error(['No post processing file found for subj' num2str(subjNum)...
                            'with weighingFactor' num2str(WeighFactor) ...
                            'Recursive Type ' RecText...
                            'Cos Type ' cosText...
                            'normalised by' normBy]);
                    end
                else
                    break;
                end
            end
            
            detailsClust = load(FinalPostProcessPath);
            voxelIndAll = detailsClust.voxelIndByCluster;
            
            for clusterNum = clustRange
                
                voxelIndByCluster1 = voxelIndAll(:,clusterNum);
                voxelIndByCluster2 = voxelIndAll(:,clusterNum+1); % this would fail if clustRange became greater than 1000
                
                ret = calcCommonClusters(voxelIndByCluster1, voxelIndByCluster2);
                
                finalret.IntersectVal{clusterNum} =  ret.IntersectVal;
                finalret.IntersectInd{clusterNum,1} = ret.IntersectInd1;
                finalret.IntersectInd{clusterNum,2} =  ret.IntersectInd2;
                finalret.IntersectCount{clusterNum} = ret.IntersectCount;
                finalret.IntersectRatio{clusterNum,1} = ret.IntersectRatio1;
                finalret.IntersectRatio{clusterNum,2} = ret.IntersectRatio2;
                finalret.MaxIntRatio{clusterNum,1} = ret.MaxIntRatio1;
                finalret.MaxIntRatio{clusterNum,2} = ret.MaxIntRatio2;
                finalret.IndMaxIntRatio{clusterNum,1} = ret.IndMaxIntRatio1;
                finalret.IndMaxIntRatio{clusterNum,2} = ret.IndMaxIntRatio2;
                
            end
        end
    end
elseif CosTypeCompare
    
    prompt = {'Enter Weighing Factor range:',...
        'Enter subject range:',...
        'Enter clustering range:',...
        'Enter Decay Parameter:',...
        'Calculate for Recursive Split("true" or "false")?',...
        'Enter individual normalisation factor("sum" or "mean")',...
        'Enter whole normalisation("WholeMax")',...
        'Enter Thresholding Factor Search Range'...
        'Enter PostProcessing Path:'};
    dlg_title = 'FS Path';
    num_lines = 1;
    defaultans = {'/net/store/nbp/projects/phasesim/workdir/Arushi/20160606_voxelByIndFS/',...
        'false',...
        '(0:0.1:1)',...
        '[1]',...
        '(66:100)',...
        '-1', ...
        'false',...
        'true',...
        'sum',...
        'WholeMax',...
        '[100,10,5,1]',...
        '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    FSPath        = answer{1};
    OverlapFlag   = str2num(answer{2});
    WeighFacRange = str2num(answer{3});
    subjRange     = str2num(answer{4});
    clustRange    = str2num(answer{5});
    decayParam    = answer{6};
    useCosineSim  = str2num(answer{7});
    recursiveSplit = str2num(answer{8});
    normBy         = answer{9};
    WholeNormText  = answer{10};
    threshRange    = str2num(answer{11});
    PostProcessPath = answer{12};
    
    if recursiveSplit
        RecText = 'Rec';
    else
        RecText = 'NonRec';
    end
    
    if useCosineSim
        cosText = 'cos';
    else
        cosText = 'conn';
    end
    
    for subjNum = subjRange
        if subjNum >= 10
            caNum = [num2str(subjNum)];
        else
            caNum = ['0' num2str(subjNum)];
        end
        for WeighFactor = WeighFacRange
            
            threshCount = 1;
            while 1
                
                FinalPostProcessPath = [PostProcessPath RecText '/decay' decayParam 'weigh'...
                    num2str(WeighFactor) '/normby' normBy 'thresh' ...
                    num2str(threshRange(threshCount)) '/' cosText '/Subj' num2str(subjNum) '/detailsClust.mat'];
                
                if ~exist(FinalPostProcessPath, 'file')
                    threshCount = threshCount + 1;
                    if threshCount > length(threshRange)
                        error(['No post processing file found for subj' num2str(subjNum)...
                            'with weighingFactor' num2str(WeighFactor) ...
                            'Recursive Type ' RecText...
                            'Cos Type ' cosText...
                            'normalised by' normBy]);
                    end
                else
                    break;
                end
            end
            
            detailsClust = load(FinalPostProcessPath);
            voxelIndAll = detailsClust.voxelIndByCluster;
            
            for clusterNum = clustRange
                
                voxelIndByCluster1 = voxelIndAll(:,clusterNum);
                voxelIndByCluster2 = voxelByInd;
                
                ret = calcCommonClusters(voxelIndByCluster1, voxelIndByCluster2);
                
                finalret.IntersectVal{clusterNum} =  ret.IntersectVal;
                finalret.IntersectInd{clusterNum,1} = ret.IntersectInd1;
                finalret.IntersectInd{clusterNum,2} =  ret.IntersectInd2;
                finalret.IntersectCount{clusterNum} = ret.IntersectCount;
                finalret.IntersectRatio{clusterNum,1} = ret.IntersectRatio1;
                finalret.IntersectRatio{clusterNum,2} = ret.IntersectRatio2;
                finalret.MaxIntRatio{clusterNum,1} = ret.MaxIntRatio1;
                finalret.MaxIntRatio{clusterNum,2} = ret.MaxIntRatio2;
                finalret.IndMaxIntRatio{clusterNum,1} = ret.IndMaxIntRatio1;
                finalret.IndMaxIntRatio{clusterNum,2} = ret.IndMaxIntRatio2;
                
            end
        end
    end
    
end







