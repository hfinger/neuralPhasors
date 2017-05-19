% 20160504 - Changed to a script to make it easier instead of having too
%            many input parameters
%          - added inputdlg

prompt = {'Enter Weighing Factor range:',...
    'Enter subject range:',...
    'Enter clustering range:',...
    'Enter Decay Parameter:',...
    'Enter Pre-Graclus Thresholding Factor:',...
    'Do you want to use Cosine Similarity("true" or "false")?',...
    'Do you want to split Recursively("true" or "false")?',...
    'Enter individual normalisation factor("sum" or "mean")',...
    'Enter whole normalisation("WholeMax")',...
    'Enter CompSimPath:',...
    'Enter GraclusPath:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'(0:0.1:1)',...
    '[1]',...
    '(2:1000)',...
    '-1', ...
    '1',...
    'true',...
    'true',...
    'sum',...
    'WholeMax',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_GraclusCut/'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

% %comment for inputdlg
% answer = defaultans;
% answer{2} = '1';
% answer{7} = 'false';

if ~isempty(answer)

WeighFacRange = str2num(answer{1});
subjRange     = str2num(answer{2});
clustRange    = str2num(answer{3});
decayParam    = str2num(answer{4});
threshFactor  = str2num(answer{5});
useCosineSim  = str2num(answer{6});
recursiveSplit= str2num(answer{7});
normBy        = answer{8};
WholeNormText = answer{9};
CompSimPath   = answer{10};
GraclusPath   = answer{11};

if recursiveSplit
    clustRange = clustRange(end);
end
if useCosineSim
    cosText = 'cos';
else
    cosText = 'conn';
end
if recursiveSplit
    splitType = 'Rec';
else
    splitType = 'NonRec';
end
for subjNum = subjRange
    for WeighingFactor = WeighFacRange
        for clusterCount = clustRange
            
            InputPath = [CompSimPath  cosText '/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor)];
            OutputPath = [GraclusPath  cosText '/' splitType '/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor) '/' ...
                num2str(subjNum) '/normby' normBy 'thresh' num2str(threshFactor) '/' ];
            disp( [cosText 'WEIGHINGFACTOR:' num2str(WeighingFactor) 'Thresh' num2str(threshFactor)]);
            
            compSimMat = load([InputPath '/compSimMat' WholeNormText '/' 'compSimMat' WholeNormText  normBy num2str(subjNum)]);
            
            compSimMat = compSimMat.compSimMat;
%             compSimMat = compSimMat*threshFactor;
                      disp(['do clustering subj:' num2str(subjNum) 'clustercount:' num2str(clusterCount) ])
                        
                        [clusterIdPerVoxel, clusterIdPerVoxelCurrent, largestClusterId, cutValue]...
                            = applyClustering( compSimMat, clusterCount, recursiveSplit );
                                         
        
            stepWiseClusters = clusterIdPerVoxel;
            allClusters = clusterIdPerVoxelCurrent;
            
            
            if ~exist(OutputPath, 'dir')
                mkdir(OutputPath);
            end
            
            
            save([OutputPath '/graclusResultClust'  num2str(clusterCount)],...
                'stepWiseClusters', 'allClusters', 'largestClusterId', 'cutValue');
            
        end
    end
end
end 