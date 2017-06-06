prompt = {'Enter Weighing Factor range:',...
          'Enter subject range:',...
          'Enter clustering range:',...
          'Enter Decay Parameter:',...
          'Do you want to use Cosine Similarity("true" or "false")?',...
          'Enter individual normalisation factor("sum" or "mean")',...
          'Enter whole normalisation("WholeMax")',...
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
    '/net/store/nbp/projects/phasesim/workdir/Arushi/NewData/20160709_GraclusCut/',...
    '[100,10,5,1]'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
WeighFacRange = str2num(answer{1});
subjRange     = str2num(answer{2});
clustRange    = str2num(answer{3});
decayParam    = str2num(answer{4});
useCosineSim  = str2num(answer{5});
normBy         = answer{6};
WholeNormText  = answer{7};
GraclusPath    = answer{8};
threshRange    = str2num(answer{9});

if useCosineSim
    cosText = 'cos';
else
    cosText = 'conn';
 end

    splitType = 'NonRec';


for subjNum = subjRange
    for WeighingFactor = WeighFacRange
        
         InputGraclusPath = [GraclusPath cosText '/' splitType '/decay'...
                             num2str(decayParam) 'weigh' num2str(WeighingFactor)...
                             '/' num2str(subjNum) '/normby' normBy 'thresh1'];
   
    clusterPath = [InputGraclusPath '/graclusResultClust'];
    
   
     cutValue = zeros(clustRange(end),1);
        
   for clustNum = clustRange
       
       finalClusterPath = [clusterPath num2str(clustNum) '.mat'];
       cluster = load(finalClusterPath);     
       
       if clustNum == clustRange(1)
           stepWiseClusters = zeros(size(cluster.stepWiseClusters,1),1);
       end     
       stepWiseClusters = horzcat(stepWiseClusters,cluster.stepWiseClusters(:));
       
       cutValue(clustNum) = cluster.cutValue;
   end
   disp(['subj' num2str(subjNum) 'weigh' num2str(WeighingFactor)]);
   save([clusterPath 'collected'], 'stepWiseClusters', 'cutValue');
    end
end
end
       
       
       
       

