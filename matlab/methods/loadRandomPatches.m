function [ patches ] = loadRandomPatches( patchDim1, patchDim2, featureDim, numPatches, pathlist, filelist, BUFF, applyConn )
%LOADRANDOMPATCHES Summary of this function goes here
%   Detailed explanation goes here

disp('initialize variable for patches...')
patches=zeros(patchDim1, patchDim2, featureDim, numPatches);

% Determine how many patches to take from each activity file
numsamplesFromEachImg = ones(size(filelist))*floor(numPatches/numel(filelist));
numremaining = numPatches-sum(numsamplesFromEachImg);
tmp = randperm(numel(filelist));
numsamplesFromEachImg(tmp(1:numremaining)) = numsamplesFromEachImg(tmp(1:numremaining)) + 1;

% extract random patches to compute whitening matrix:
totalsamples = 0;
disp('load images and fill the patches...')
for fileid=1:numel(filelist)
%   disp(num2str(fileid))
  act = load(fullfile(pathlist{fileid},filelist{fileid}));
%   disp(fullfile(pathlist{fileid},filelist{fileid}));
  act = act.act;
  
  if nargin>7 && ~isempty(applyConn)
    act = ApplyWeights.applyWeights( act, [], applyConn);
  end
  
  % Extract patches at random from this image to make data vector X
  for k=1:numsamplesFromEachImg(fileid)
    totalsamples = totalsamples + 1;
    if iscell(act)
      cellid = ceil(numel(act) * k / numsamplesFromEachImg(fileid));
      patches(:,:,:,totalsamples) = selectRandPatch(act{cellid}, patchDim1, patchDim2, BUFF);
    else
      patches(:,:,:,totalsamples) = selectRandPatch(act, patchDim1, patchDim2, BUFF);
    end
  end
end

end

function patch = selectRandPatch(act, patchDim1, patchDim2, BUFF)
dim1start=1+BUFF+floor((size(act,1)-patchDim1-2*BUFF)*rand);
dim2start=1+BUFF+floor((size(act,2)-patchDim2-2*BUFF)*rand);
patch = act(...
  dim1start:dim1start+patchDim1-1,...
  dim2start:dim2start+patchDim2-1,...
  :);
end

