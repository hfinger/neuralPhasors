loadpath = 'noRenorm_noConstrain_noReg_moreNoiseInterp_runLonger';

taskIds=1:2;
plot2D_mds = false;

allSC = cell(1,length(taskIds));
allIters = cell(1,length(taskIds));
allNoise = cell(1,length(taskIds));
allRandIds = cell(1,length(taskIds));
allDev = cell(1,length(taskIds));
for i=1:length(taskIds)
  
  %% load data:
  tmp = load(fullfile(loadpath,['taskId' num2str(taskIds(i)) '.mat']));
  logLearnSC = cellfun(@(x) cell2mat(permute(x,[2 3 1])), tmp.logLearnSC, 'UniformOutput', false);
  
  iters = cell(size(logLearnSC));
  noise = cell(size(logLearnSC));
  dev = cell(size(logLearnSC));
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
    noise{m} = tmp.noiseAmount(m) * ones(size(iters{m}));
    
    %% save deviance at the corresponding positions:
    % add 1 because we gathered the dev within the next run.
    dev{m} = tmp.dev{m}(min(tmp.p.saveAtIters+1,length(tmp.dev{m})));
    
  end
  
  allSC{i} = cat(3,logLearnSC{:});
  allIters{i} = cat(2,iters{:});
  allNoise{i} = cat(2,noise{:});
  allRandIds{i} = i * ones(size(allNoise{i}));
  allDev{i} = dev;
  
end

allSC = cat(3,allSC{:});
allIters = cat(2,allIters{:});
allNoise = cat(2,allNoise{:});
allRandIds = cat(2,allRandIds{:});


%%
vecSC = reshape(allSC,[66*66 size(allSC,3)]);
offDiagIds = find(~eye(66,66));
vecSC = vecSC(offDiagIds,:);

%% calc similarity matrix:

finalLearnSC = zeros(length(offDiagIds),length(taskIds),length(tmp.finalLearnSC));
for i=1:length(taskIds)
  tmp = load(fullfile(loadpath,['taskId' num2str(taskIds(i)) '.mat']));
  for m=1:length(tmp.finalLearnSC)
    finalLearnSC(:,i,m) = tmp.finalLearnSC{m}(offDiagIds);
  end
end

figure(1)
distances = pdist(reshape(finalLearnSC,[size(finalLearnSC,1) size(finalLearnSC,2)*size(finalLearnSC,3)])','correlation');
imagesc(squareform(distances))
title('correlation dist between learned SC')

figure(2)
distances = pdist(reshape(finalLearnSC,[size(finalLearnSC,1) size(finalLearnSC,2)*size(finalLearnSC,3)])','euclidean');
imagesc(squareform(distances))
title('euclidean dist between learned SC')

figure(3)
distances = pdist(reshape(finalLearnSC,[size(finalLearnSC,1) size(finalLearnSC,2)*size(finalLearnSC,3)])','cosine');
imagesc(squareform(distances))
title('cosine dist between learned SC')


%% find all empSC and remove them so that we dont have duplicates:
% empSCids = find(sum(abs(bsxfun(@minus, vecSC(:,1), vecSC )), 1)==0);
% numRemoved = length(empSCids)-1;
% 
% tmpRemovedVecSC = vecSC(:,empSCids(2:end));
% tmpRemovedSC = allSC(:,:,empSCids(2:end));
% tmpRemovedIters = allIters(empSCids(2:end));
% tmpRemovedNoise = allNoise(empSCids(2:end));
% tmpRemovedRandIds = allRandIds(empSCids(2:end));
% 
% vecSC(:,empSCids(2:end)) = [];
% allSC(:,:,empSCids(2:end)) = [];
% allIters(empSCids(2:end)) = [];
% allNoise(empSCids(2:end)) = [];
% allRandIds(empSCids(2:end)) = [];

%% 
distances = pdist(vecSC','correlation');
zeroDistIds = find(distances==0);

%%
if plot2D_mds
  Y2 = mdscale(distances,2,'Criterion','sstress');
end
Y3 = mdscale(distances,3,'Criterion','sstress');

%% append removed items:
% Y2(end+1:end+numRemoved,:) = repmat(Y2(1,:),[numRemoved, 1]);
% Y3(end+1:end+numRemoved,:) = repmat(Y3(1,:),[numRemoved, 1]);
% vecSC(:,end+1:end+numRemoved) = tmpRemovedVecSC;
% allSC(:,:,end+1:end+numRemoved) = tmpRemovedSC;
% allIters(end+1:end+numRemoved) = tmpRemovedIters;
% allNoise(end+1:end+numRemoved) = tmpRemovedNoise;
% allRandIds(end+1:end+numRemoved) = tmpRemovedRandIds;


%%
if plot2D_mds
  figure(4)
  plot(Y2(:,1),Y2(:,2),'.')
  title('2D MDS')
end

figure(5)
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'.')
title('3D MDS')

%%
if plot2D_mds
  figure(6)
  clf;
  hold on;
  for i=1:length(taskIds)
    for m=1:length(logLearnSC)
      SCids = find(allRandIds==taskIds(i) & allNoise==tmp.noiseAmount(m));
      [~, I] = sort(allIters(SCids));
      SCids = SCids(I);
      plot(Y2(SCids,1),Y2(SCids,2),'r')
    end
    SCids = find(allRandIds==taskIds(i) & allIters==0);
    [~, I] = sort(allNoise(SCids));
    SCids = SCids(I);
    plot(Y2(SCids,1),Y2(SCids,2))
  end
  plot(Y2(1,1),Y2(1,2),'bo')
  title('2D MDS')
end

%% plot 3d mds
figure(7)
clf;
hold on;
for i=1:length(taskIds)
  for m=1:length(logLearnSC)
    SCids = find(allRandIds==taskIds(i) & allNoise==tmp.noiseAmount(m));
    [~, I] = sort(allIters(SCids));
    SCids = SCids(I);
    plot3(Y3(SCids,1),Y3(SCids,2),Y3(SCids,3),'.-r')
  end
  SCids = find(allRandIds==taskIds(i) & allIters==0);
  [~, I] = sort(allNoise(SCids));
  SCids = SCids(I);
  plot3(Y3(SCids,1),Y3(SCids,2),Y3(SCids,3))
end
plot3(Y3(1,1),Y3(1,2),Y3(1,3),'bo')
title('3D MDS')
axis equal;


%% plot 3d mds with color as loss:
figure(7)
clf;
colormap(jet)
hold on;
for i=1:length(taskIds)
  for m=1:length(logLearnSC)
    SCids = find(allRandIds==taskIds(i) & allNoise==tmp.noiseAmount(m));
    [~, I] = sort(allIters(SCids));
    SCids = SCids(I);
    
    %% search color:
    my_dev = allDev{i}{m};
    
    surface([Y3(SCids,1)';Y3(SCids,1)'],[Y3(SCids,2)';Y3(SCids,2)'],[Y3(SCids,3)';Y3(SCids,3)'],[my_dev';my_dev'],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
    
    %plot3(Y3(SCids,1),Y3(SCids,2),Y3(SCids,3),'.-r')
  end
  SCids = find(allRandIds==taskIds(i) & allIters==0);
  [~, I] = sort(allNoise(SCids));
  SCids = SCids(I);
  plot3(Y3(SCids,1),Y3(SCids,2),Y3(SCids,3),'linew',2,'Color','k')
end
plot3(Y3(1,1),Y3(1,2),Y3(1,3),'bo')
title('3D MDS')
colorbar
set(gca,'clim',[0 1])
axis equal;
