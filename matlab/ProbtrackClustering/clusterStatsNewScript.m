% clusterStatsNewScript

prompt = {'Enter Weighing Factor range:',...
    'Enter subject range:',...
    'Enter clustering range:',...
    'Enter Decay Parameter:',...
    'Do you want to use Cosine Similarity("true" or "false")?',...
    'Calculate for Recursive Split("true" or "false")?',...
    'Enter individual normalisation factor("sum" or "mean")',...
    'Enter whole normalisation("WholeMax")',...
    'Enter PostProcessing Path',...
    'Enter Output Path',...
    'Enter Thresholding Factor Search Range'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'(0:0.1:1)',...
    '[1]',...
    '(2:1000)',...
    '-1', ...
    'false',...
    'true',...
    'sum',...
    'WholeMax',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/'...
    '[100,10,5,1]'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
    WeighFacRange = str2num(answer{1});
    subjRange     = str2num(answer{2});
    clustRange    = str2num(answer{3});
    decayParam    = str2num(answer{4});
    useCosineSim  = str2num(answer{5});
    recursiveSplit = str2num(answer{6});
    normBy         = answer{7};
    WholeNormText  = answer{8};
    GraclusPath    = answer{9};
    PostProcessPath = answer{10};
    threshRange    = str2num(answer{11});
    
    if useCosineSim
        cosText = 'cos';
        SFinalType = Stype;
    else
        cosText = 'conn';
        SFinalType = 'SSym';
    end
    if recursiveSplit
        splitType = 'Rec';
    else
        splitType = 'NonRec';
    end
    
    
    for subjNum = subjRange
        for WeighingFactor = WeighFacRange
            
            disp(['subj' num2str(subjNum) 'weigh' num2str(WeighingFactor)]);
            for threshFactor = threshRange
                if exist([PostProcessPath splitType '/decay' num2str(decayParam)...
                'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh'...
                num2str(threshFactor) '/' cosText  '/' 'Subj'...
                num2str(subjNum) '/'], 'dir')
                    threshFound = 1;
                    break;
                end
            end
            
            if ~threshFound
                error('No cluster matrix found for subj %i for weighingFactor %i',subjNum, WeighingFactor);
            end
            clusterPath = [PostProcessPath splitType '/decay' num2str(decayParam)...
                'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh'...
                num2str(threshFactor) '/' cosText  '/' 'Subj'...
                num2str(subjNum) '/'];
           
            %      clusterPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/Rec/ClustOrgW0/detailsClust'];
            %     OutputPath = ['/net/store/nbp/projects/phasesim/workdir/Arushi/20160422_ClusteringPostprocessing/Subj' num2str(subjNum) '/Rec/ClustOrgW0/StatandMet/'];
          
            clusterStatsNew(clustRange, clusterPath, clusterPath);
            
            
        end
    end
end
