function [ H ] = linearDecoderNghWinMat( topoDimNgh,topoDimSize,topoDimPeriodicBoundary, topoNghWin )
%LINEARDECODERNGHMAT Summary of this function goes here
%   Detailed explanation goes here

% Example:
% topoDimNgh = [3 3];
% topoDimSize = [20 20];
% topoDimPeriodicBoundary = [true true];


  
if length(topoDimNgh) == 1;
  topoDimNgh = repmat(topoDimNgh,[1 length(topoDimSize)]);
end
if length(topoDimPeriodicBoundary) == 1;
  topoDimPeriodicBoundary = repmat(topoDimPeriodicBoundary,[1 length(topoDimSize)]);
end

numNeurons = prod(topoDimSize);
H=zeros([numNeurons numNeurons]);
coord = cell(size(topoDimSize));

if length(topoDimSize)==1
  nghWin = ones(2*topoDimSize+1,1);
else
  nghWin = ones(2*topoDimSize+1);
end

for dimId=1:length(topoDimSize)
%   w = feval(topoNghWin,2*topoDimSize(dimId)+1,topoDimNgh(dimId));
  w = feval(topoNghWin,2*topoDimSize(dimId)+1);
  topoDimSizeTmp=2*topoDimSize+1;
  topoDimSizeTmp(dimId)=1;
  tmp=ones(size(topoDimSize));
  tmp(dimId) = numel(w);
  if length(tmp)==1
    tmp=[tmp(1) 1];
    topoDimSizeTmp = [1 1];
  end
  nghWin = nghWin .* repmat(reshape(w,tmp),topoDimSizeTmp);
end

idxInDims = ones(size(topoDimSize));
idxInDims(1)=0;
for i=1:numNeurons
  %% find next sub dims (like ind2sub):
  idxInDims(1)=idxInDims(1)+1;
  for dimId=1:length(topoDimSize)
    if idxInDims(dimId)>topoDimSize(dimId)
      idxInDims(dimId)=1;
      idxInDims(dimId+1)=idxInDims(dimId+1)+1;
    end
  end
  
  %% find neighboors:
  for dimId=1:length(topoDimSize)
    tmpCoord = (topoDimSize(dimId)-idxInDims(dimId)+2):(topoDimSize(dimId)-idxInDims(dimId)+topoDimSize(dimId)+1);
    if topoDimPeriodicBoundary(dimId)
      tmpCoord(tmpCoord<=topoDimSize(dimId)/2)=tmpCoord(tmpCoord<=topoDimSize(dimId)/2)+topoDimSize(dimId);
      tmpCoord(tmpCoord>3*topoDimSize(dimId)/2)=tmpCoord(tmpCoord>3*topoDimSize(dimId)/2)-topoDimSize(dimId);
    end
    coord{dimId} = tmpCoord;
  end
  H( i, : ) = reshape(nghWin(coord{:}),[1 size(H,2)]);
end

H = reshape(H,[numNeurons numNeurons]);
end

