function img = plotImages( imgCell, savefile, padLength )
%PLOTIMAGES Summary of this function goes here
%   Detailed explanation goes here

if ~iscell(imgCell)
  % imgCell has to be a 5 dimensional array with [x,y,rgb,patchIdX,patchIdY]
  % convert dimensions 4 and 5 to cell:
  imgCell=squeeze(num2cell(imgCell,[1 2 3]));
end

imgCell = cellfun(@(x) cat(2,x,255*ones(size(x,1),padLength,3)), imgCell,'UniformOutput',false);
imgCell = cellfun(@(x) cat(1,x,255*ones(padLength,size(x,2),3)), imgCell,'UniformOutput',false);

img=cell2mat(imgCell);
img=img(:,1:end-padLength,:);
img=img(1:end-padLength,:,:);

if ~isempty(savefile)
  imwrite(img,savefile,'png')
end

end

