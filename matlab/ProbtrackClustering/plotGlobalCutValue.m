%GlobalCutValCalcscript

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
defaultans = {'(0:0.5:1)',...
  '[1]',...
  '(2:1000)',...
  '-1', ...
  'false',...
  'sum',...
  'WholeMax',...
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
  GraclusPath    = answer{8};
  threshRange    = str2num(answer{9});
  
  if useCosineSim
    cosText = 'cos';
  else
    cosText = 'conn';
  end
  
  FSCutVal = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20160706_FScalc/' cosText '/FSCutVal.mat']);
  

  
  
  for subjNum = subjRange
    
    figure(1);
    clf;
    hold on;
    set(gca,'FontSize',16)
    set(gca,'LineWidth',2)
  
    k=1;
    for WeighingFactor = WeighFacRange(end:-1:1)
      disp(['subj' num2str(subjNum) 'weigh' num2str(WeighingFactor)]);
      
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
      
      data = load([GraclusPath 'globalCutVal/Rec' cosText '/subj' num2str(subjNum) '/decay'...
        num2str(decayParam) 'weigh' num2str(WeighingFactor) '/normby' normBy 'thresh' num2str(threshFactor) '/GlobalCutVal.mat']);
      
      h(k) = plot( data.cumulGlobalCutVal ./ (1:1000), 'LineWidth', 2 );
      
      ax = gca;
      ax.ColorOrderIndex = ax.ColorOrderIndex-1;
      plot(66, FSCutVal.cumulGlobalCutVal(1+round(WeighingFactor*10))/66, 'x', 'MarkerSize', 10, 'LineWidth', 2);
      
      k=k+1;
      
    end
    
    h(end+1) = plot(-Inf, -Inf, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    
    legend(h,{'Local','DTI+Local','DTI','Freesurfer'}, 'Location', 'NorthWest')
    ylabel('Average Cut-Value')
    xlabel('Number of Clusters')
  
  end
end