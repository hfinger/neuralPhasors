function img = plotColorFeatures( img, scalePerPatch, output, divMax, padLength )
%PLOTCOLORFEATURES Summary of this function goes here
%   Detailed explanation goes here

%img should be dimx,dimy,3,numimagesX,numImagesY

if nargin<4
  divMax = false;
end
if nargin<5
  padLength = 1;
end

% imgCell=squeeze(num2cell(img,[1 2 3]));
imgCell=permute(num2cell(img,[1 2 3]),[4 5 1 2 3]);
if divMax
  if scalePerPatch
    imgCell = cellfun(@(x) ((x/max(abs(x(:))))+1)/2, imgCell,'UniformOutput',false);
  else
    imgCell = cellfun(@(x) ((x/max(abs(img(:))))+1)/2, imgCell,'UniformOutput',false);
  end
else
  if scalePerPatch
    imgCell = cellfun(@(x) (x-min(x(:)))/(max(x(:))-min(x(:))), imgCell,'UniformOutput',false);
  else
    imgCell = cellfun(@(x) (x-min(img(:)))/(max(img(:))-min(img(:))), imgCell,'UniformOutput',false);
  end
end
img = plotImages( imgCell', output, padLength );

end

