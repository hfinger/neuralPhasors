function [img] = plotFeatures( features, savefile, strColormap, padLength, cmin, cmax, climPi, genPdf, noPlot, invertColormap )
%PLOTFEATURES Summary of this function goes here
%   Detailed explanation goes here

if nargin<4 || isempty(padLength)
  padLength=0;
end
if nargin<3 || isempty(strColormap)
  strColormap='jet';
end
if nargin<2 || isempty(savefile)
  savefile = [];
end
if nargin<8 || isempty(genPdf)
  genPdf = false;
end
if nargin<9 || isempty(noPlot)
  noPlot = false;
end
if nargin<10 || isempty(invertColormap)
  invertColormap = false;
end

if iscell(features)
  features = cat(4,features{:});
  %   for k=1:length(features)
  %     dim=size(features{k});
  %     features{k} = reshape(features{k}, dim(1), dim(2), prod(dim(3:end)) );
  %   end
  %   features = cat(3,features{:});
end

if nargin<5 || isempty(cmin)
  a=min(features(:));
else
  a=cmin;
  features(features<a)=a;
end

if nargin<6 || isempty(cmax)
  b=max(features(:));
else
  b=cmax;
  features(features>b)=b;
end

if nargin<7 || isempty(climPi)
  climPi=[];
end

if invertColormap
  cmap = colormap(flipud(colormap(strColormap)));
else
  cmap = colormap(strColormap);
end

% features = uint8( (features-a)*256/(b-a) );
features = (features-a)/(b-a);

imgCell = num2cell(features, [1 2]);
imgCell = permute(imgCell,[3 4 1 2]);
imgCell = cellfun(@(x) grs2rgb(x,cmap,1), imgCell,'UniformOutput',false);
imgCell = cellfun(@(x) cat(2,x,ones(size(x,1),padLength,3)), imgCell,'UniformOutput',false);
imgCell = cellfun(@(x) cat(1,x,ones(padLength,size(x,2),3)), imgCell,'UniformOutput',false);

% img=cat(2,imgCell{:});
img=cell2mat(imgCell);
img=img(:,1:end-padLength,:);
img=img(1:end-padLength,:,:);
%imwrite(img,savefile,'png');

%% export colorbar:
% hFig=figure('Visible','off');
% colormap(strColormap);
% caxis([a b]);
% hC=colorbar;
% rect=getframe(hC);
% close(hFig);
% tmpImg=frame2im(rect);
% imwrite(tmpImg,[savefile 'ColorbarMin' num2str(a) 'Max' num2str(b) '.png'],'png');

if ~noPlot
  %% export as tikz:
  if a==b
    a=a-0.1;
    b=b+0.1;
  end
  
  imagesc(img,[a b])
  colormap(cmap);
  hC=colorbar;
  axis image;
  if ~isempty(savefile)
    exportTikZImage( savefile, climPi, genPdf );
  end
end
end

