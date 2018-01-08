taskIds=1:10;

allSC = cell(1,length(taskIds));
allIters = cell(1,length(taskIds));
allNoise = cell(1,length(taskIds));
allRandIds = cell(1,length(taskIds));
for i=1:length(taskIds)
  
  %% load data:
  tmp = load(['taskId' num2str(taskIds(i)) '.mat']);
  logLearnSC = cellfun(@(x) cell2mat(permute(x,[2 3 1])), tmp.logLearnSC, 'UniformOutput', false);
  
  iters = cell(size(logLearnSC));
  noise = cell(size(logLearnSC));
  for m=1:length(logLearnSC)
    %%
    numSaved = size(logLearnSC{m},3);
    iters{m} = tmp.p.saveAtIters(1:numSaved);
    
    %% append last iter
    lastIter = find(tmp.dev{m}==0,1,'first');
    if ~isempty(lastIter)
      logLearnSC{m}(:,:,end+1) = tmp.finalLearnSC{m};
      iters{m}(end+1) = lastIter;
    end
    noise{m} = tmp.interpPoints(m) * ones(size(iters{m}));
  
  end
  
  allSC{i} = cat(3,logLearnSC{:});
  allIters{i} = cat(2,iters{:});
  allNoise{i} = cat(2,noise{:});
  allRandIds{i} = i * ones(size(allNoise{i}));
  
end

allSC = cat(3,allSC{:});
allIters = cat(2,allIters{:});
allNoise = cat(2,allNoise{:});
allRandIds = cat(2,allRandIds{:});

%% keep only full noise and one with no noise:
keepIds = [find(allNoise==1) find(allNoise==0 & allRandIds==1)];

allSC = allSC(:,:,keepIds);
allIters = allIters(keepIds);
allNoise = allNoise(keepIds);
allRandIds = allRandIds(keepIds);

%%
vecSC = reshape(allSC,[66*66 size(allSC,3)]);
offDiagIds = find(~eye(66,66));
vecSC = vecSC(offDiagIds,:);

%% find all empSC and remove them so that we dont have duplicates:
empSCids = find(sum(abs(bsxfun(@minus, vecSC(:,1), vecSC )), 1)==0);
numRemoved = length(empSCids)-1;

tmpRemovedVecSC = vecSC(:,empSCids(2:end));
tmpRemovedSC = allSC(:,:,empSCids(2:end));
tmpRemovedIters = allIters(empSCids(2:end));
tmpRemovedNoise = allNoise(empSCids(2:end));
tmpRemovedRandIds = allRandIds(empSCids(2:end));

vecSC(:,empSCids(2:end)) = [];
allSC(:,:,empSCids(2:end)) = [];
allIters(empSCids(2:end)) = [];
allNoise(empSCids(2:end)) = [];
allRandIds(empSCids(2:end)) = [];

%% 
distances = pdist(vecSC','correlation');
zeroDistIds = find(distances==0);

%%
Y2 = mdscale(distances,2,'Criterion','sstress');
Y3 = mdscale(distances,3,'Criterion','sstress');

%% append removed items:
Y2(end+1:end+numRemoved,:) = repmat(Y2(1,:),[numRemoved, 1]);
Y3(end+1:end+numRemoved,:) = repmat(Y3(1,:),[numRemoved, 1]);
vecSC(:,end+1:end+numRemoved) = tmpRemovedVecSC;
allSC(:,:,end+1:end+numRemoved) = tmpRemovedSC;
allIters(end+1:end+numRemoved) = tmpRemovedIters;
allNoise(end+1:end+numRemoved) = tmpRemovedNoise;
allRandIds(end+1:end+numRemoved) = tmpRemovedRandIds;


%%

figure(1)
plot(Y2(:,1),Y2(:,2),'.')

figure(2)
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'.')

%%
figure(3)
clf;
hold on;
for i=1:length(taskIds)
  for m=1:length(logLearnSC)
    SCids = find(allRandIds==taskIds(i) & allNoise==tmp.interpPoints(m));
    [~, I] = sort(allIters(SCids));
    SCids = SCids(I);
    plot(Y2(SCids,1),Y2(SCids,2),'r')
  end
  SCids = find(allRandIds==taskIds(i) & allIters==0);
  [~, I] = sort(allNoise(SCids));
  SCids = SCids(I);
  plot(Y2(SCids,1),Y2(SCids,2),'b')
end
plot(Y2(1,1),Y2(1,2),'bo')
