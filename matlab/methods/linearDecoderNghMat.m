function [ H ] = linearDecoderNghMat( topoDimNgh,topoDimSize,topoDimPeriodicBoundary )
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

H=zeros([numNeurons topoDimSize]);
coord = cell(size(topoDimSize));

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
  coord{1} = i;
  for dimId=1:length(topoDimSize)
    tmpCoord = (idxInDims(dimId)-floor(topoDimNgh(dimId)/2)):(idxInDims(dimId)+ceil(topoDimNgh(dimId)/2)-1);
    if topoDimPeriodicBoundary(dimId)
      tmpCoord(tmpCoord<1)=tmpCoord(tmpCoord<1)+topoDimSize(dimId);
      tmpCoord(tmpCoord>topoDimSize(dimId))=tmpCoord(tmpCoord>topoDimSize(dimId))-topoDimSize(dimId);
    else
      tmpCoord(tmpCoord<1)=[];
      tmpCoord(tmpCoord>topoDimSize(dimId))=[];
    end
    coord{dimId+1} = tmpCoord;
  end
  H( coord{:} ) = 1;
end

H = reshape(H,[numNeurons numNeurons]);
end

