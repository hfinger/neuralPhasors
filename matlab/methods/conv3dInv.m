function [ W ] = conv3dInv( z2,data )
%CONV3DINV Pairwise combination between z2 and data. 
% Assuming that z2(j,k,f2) = W(x,y,f1,f2) * data(x+j,y+k,f1)

if size(z2,1)==1 && size(z2,2)==1 && size(data,1)==1 && size(data,2)==1
  W = reshape( reshape(data,[size(data,3) 1]) * reshape(z2,[1 size(z2,3)]), [1 1 size(data,3) size(z2,3)]);
else
  W = zeros([size(data,1)-size(z2,1)+1 size(data,2)-size(z2,2)+1 size(data,3) size(z2,3)]);
  for f1=1:size(data,3)
    blockImg = im2col(data(:,:,f1),[size(z2,1) size(z2,2)]);
    blockFilters = reshape(z2(:,:,:),[size(z2,1)*size(z2,2) size(z2,3)]);
    W(:,:,f1,:) = reshape(blockImg'*blockFilters,[size(W,1) size(W,2) 1 size(z2,3)]);
  end
end

end

