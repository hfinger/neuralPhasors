function kurtosis = activityKurtosis( inActFolder, inActFilenames, actHistFilename )
%ACTIVITYHIST Summary of this function goes here
%   Detailed explanation goes here
[ pathlist, filelist ] = dirrec( inActFolder , inActFilenames );

meanAct = load(actHistFilename);
meanAct = meanAct.firstMoment/meanAct.numberUnitsAdded;
numFeatures = numel(meanAct);

fourthMoment = zeros(1,numFeatures);
secondMoment = zeros(1,numFeatures);
numberUnitsAdded = 0;
for fileid=1:length(pathlist)
  act = load(fullfile(pathlist{fileid},filelist{fileid}));
  act = act.act;
  if iscell(act)
    for cellid=1:length(act)
      act{cellid} = reshape(act{cellid},[size(act{cellid},1)*size(act{cellid},2) size(act{cellid},3)]);
      act{cellid} = bsxfun(@minus,act{cellid},meanAct);
      fourthMoment = fourthMoment + sum(act{cellid}.^4,1);
      secondMoment = secondMoment + sum(act{cellid}.^2,1);
      numberUnitsAdded = numberUnitsAdded + size(act{cellid},1);
    end
  else
    act = reshape(act,[size(act,1)*size(act,2) size(act,3)]);
    act = bsxfun(@minus,act,meanAct);
    fourthMoment = fourthMoment + sum(act.^4,1);
    secondMoment = secondMoment + sum(act.^2,1);
    numberUnitsAdded = numberUnitsAdded + size(act,1);
  end
end

kurtosis = (fourthMoment/numberUnitsAdded) ./ (secondMoment/numberUnitsAdded).^2 - 3;

save([actHistFilename 'Kurtosis.mat'],'fourthMoment','secondMoment','numberUnitsAdded','kurtosis');

end

