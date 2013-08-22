function activityHist( inActFolder, inActFilenames, outSaveFile, edges, numFeatures, scaleAct, edgesScaled )
%ACTIVITYHIST Summary of this function goes here
%   Detailed explanation goes here
[ pathlist, filelist ] = dirrec( inActFolder , inActFilenames );

disp(['pwd: ' pwd]);
disp(['numfiles: ' num2str(length(filelist))])

n = zeros(length(edges),numFeatures);
nScaled = zeros(length(edgesScaled),numFeatures);
fourthMoment = zeros(1,numFeatures);
secondMoment = zeros(1,numFeatures);
firstMoment = zeros(1,numFeatures);
numberUnitsAdded = 0;
for fileid=1:length(pathlist)
  act = load(fullfile(pathlist{fileid},filelist{fileid}));
  act = act.act;
  if iscell(act)
    for cellid=1:length(act)
      act{cellid} = reshape(act{cellid},[size(act{cellid},1)*size(act{cellid},2) size(act{cellid},3)]);
      n = n + histc(act{cellid},edges);
      fourthMoment = fourthMoment + sum(act{cellid}.^4,1);
      secondMoment = secondMoment + sum(act{cellid}.^2,1);
      firstMoment = firstMoment + sum(act{cellid},1);
      numberUnitsAdded = numberUnitsAdded + size(act{cellid},1);
    end
  else
    act = reshape(act,[size(act,1)*size(act,2) size(act,3)]);
    n = n + histc(act,edges);
    fourthMoment = fourthMoment + sum(act.^4,1);
    secondMoment = secondMoment + sum(act.^2,1);
    firstMoment = firstMoment + sum(act,1);
    numberUnitsAdded = numberUnitsAdded + size(act,1);
  end
end
meanAct = firstMoment/numberUnitsAdded;

if scaleAct
  for fileid=1:length(pathlist)
    act = load(fullfile(pathlist{fileid},filelist{fileid}));
    act = act.act;
    if iscell(act)
      for cellid=1:length(act)
        act{cellid} = bsxfun(@rdivide,reshape(act{cellid},[size(act{cellid},1)*size(act{cellid},2) size(act{cellid},3)]),meanAct);
        nScaled = nScaled + histc(act{cellid},edgesScaled);
      end
    else
      act = bsxfun(@rdivide,reshape(act,[size(act,1)*size(act,2) size(act,3)]),meanAct);
      nScaled = nScaled + histc(act,edgesScaled);
    end
  end
end

save(outSaveFile,'n', 'nScaled', 'fourthMoment', 'secondMoment', 'numberUnitsAdded', 'firstMoment','edgesScaled','edges');

end

