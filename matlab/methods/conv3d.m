function z2 = conv3d(W,data,invertW)
% Inputs:
%data has dimensions x,y,fin
%W has dimensions dx,dy,fin,fout
%z2 has dimensions x,y,fout
%
% WARNING: The convolution is not flipped. So it is actually a correlation!
% so z2(j,k) = W(x,y) * data(x+j,y+k)

if nargin<3
  invertW = false;
end
if invertW
  W = permute(W(end:-1:1,end:-1:1,:,:),[1 2 4 3]); % This is actually faster than
  %W = permute(flipdim(flipdim(W,1),2),[1 2 4 3]);
end

ndx = size(W,1);
ndy = size(W,2);
nfin = size(W,3);
nfout = size(W,4);

if size(W,1)==1 && size(W,2)==1
  z2 = reshape( reshape(data,[1 size(data,3)]) * reshape(W,[nfin nfout]) , [1 1 nfout]);
else
  z2 = zeros([size(data,1)-size(W,1)+1 size(data,2)-size(W,2)+1 nfout]);
  for finId=1:size(data,3)
    blockImg = im2col(data(:,:,finId),[size(W,1) size(W,2)]);
    blockFilters = reshape(W(:,:,finId,:),[size(W,1)*size(W,2) nfout]);
    z2 = z2 + reshape(blockImg'*blockFilters,[size(data,1)-size(W,1)+1 size(data,2)-size(W,2)+1 nfout]);
  end
%   if true
%     blkSize = [size(W,1) size(W,2)];
%     W = reshape(W,[size(W,1)*size(W,2)*size(W,3) size(W,4)]);
%     %% blockproc cannot use sliding windows!
%     z2test = blockproc(data,blkSize,@(x) reshape(x.data(:)'*W,[1 1 size(W,2)]), 'PadPartialBlocks',true);
%     %% nlfilter works only for 2D and cannot process a feature dimension:
%     z2test = nlfilter(data,blkSize,@(x) reshape(x.data(:)'*W,[1 1 size(W,2)]));
%   end
end
end