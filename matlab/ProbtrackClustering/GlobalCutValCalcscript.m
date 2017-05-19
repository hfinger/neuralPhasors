%GlobalCutValCalcscript

prompt = {'Enter Weighing Factor range:',...
          'Enter subject range:',...
          'Enter clustering range:',...
          'Enter Decay Parameter:',...
          'Do you want to use Cosine Similarity("true" or "false")?',...
          'Enter individual normalisation factor("sum" or "mean")',...
          'Enter whole normalisation("WholeMax")',...
          'Enter CompSimPath:',...
          'Enter GraclusPath:',... 
          'Enter Thresholding Factor Search Range'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'(0:0.1:1)',...
    '[1]',...
    '(2:1000)',...
    '-1', ...
    'false',...
    'sum',...
    'WholeMax',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/20160418_CompSimMatCalcNewSub/',...
    '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData_copy/20160709_GraclusCut/',...
    '[1]'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
WeighFacRange = str2num(answer{1});
subjRange     = str2num(answer{2});
clustRange    = str2num(answer{3});
decayParam    = str2num(answer{4});
useCosineSim  = str2num(answer{5});
normBy         = answer{6};
WholeNormText  = answer{7};
CompSimPath    = answer{8};
GraclusPath    = answer{9};
threshRange    = str2num(answer{10});

if useCosineSim
    cosText = 'cos';
else
    cosText = 'conn';
end

for subjNum = subjRange
    for WeighingFactor = WeighFacRange
        disp(['subj' num2str(subjNum) 'weigh' num2str(WeighingFactor)]);
        
        InputCompPath = [CompSimPath  cosText '/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor) '/compSimMat' WholeNormText '/' 'compSimMat' WholeNormText  normBy num2str(subjNum)];
      InputGraclusPath = [GraclusPath  cosText '/Rec/decay'...
                num2str(decayParam) 'weigh' num2str(WeighingFactor) '/' ...
               num2str(subjNum) '/normby' normBy 'thresh'   ];
           threshFound = 0;
           for threshFactor = threshRange
               if exist([InputGraclusPath num2str(threshFactor)], 'dir')
                   threshFound = 1;
                   break;
               end
           end
           
           if ~threshFound
               error('No cluster matrix found for subj %i for weighingFactor %i',subjNum, WeighingFactor);
           end
               
                 
           
           OutputPath = [GraclusPath 'globalCutVal/Rec' cosText '/subj' num2str(subjNum) '/decay'...
               num2str(decayParam) 'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh' num2str(threshFactor) '/'];
           ClusterIdPath = [InputGraclusPath  num2str(threshFactor) '/graclusResultClust' num2str(clustRange(end))];
           
           GlobalCutValCalc( clustRange, InputCompPath, OutputPath, ClusterIdPath, threshFactor );
    end
end
           
           
end