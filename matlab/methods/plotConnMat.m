function plotConnMat( connMat, plotDir, fname, clims, upscale, onlyTriu, flipRH, doResort )
%PLOTCONNMAT Summary of this function goes here
%   Detailed explanation goes here

if nargin<5
  upscale = false;
end

if nargin<6
  onlyTriu = false;
end

if nargin<7
  flipRH = false;
end

if nargin<8
  doResort = false;
end

figh = sfigure(1);
clf;
set(figh,'Position', [100, 100, 200, 200])

numNodes = size(connMat,1);

if doResort
  paths = dataPaths();
  resortIds = load([paths.databases '/SC_Bastian/resortIds.mat']);
  resortIds = resortIds.resortIds;
  connMat = connMat(resortIds,resortIds);
end

if flipRH
  connMat(numNodes/2+1:end,:) = connMat(end:-1:numNodes/2+1,:);
  connMat(:,numNodes/2+1:end) = connMat(:,end:-1:numNodes/2+1);
end

if onlyTriu
  trilIds = find(tril(ones(66,66)));
  connMat(trilIds)=NaN; %#ok<FNDSB>
end

if upscale
  connMatHiRes = repmat(connMat,[1 1 10 10]);
  connMatHiRes = reshape(permute(connMatHiRes,[3 1 4 2]),10*size(connMat));
  imagesc(connMatHiRes)
else
  imagesc(connMat)
end
axis image;

% title(['average subject SC'])
% colormap hot;

c = colorbar;
if nargin>3 && ~isempty(clims)
  set(gca,'clim',clims);
end

% set(c,'YTick',get(gca,'clim'));

x1=get(gca,'position');
x=get(c,'Position');
x(3)=0.03;
set(c,'Position',x)
x1(1)=0.04+x1(1);
set(gca,'position',x1)

ticPos = [16 49];
if upscale
  ticPos = ticPos * 10;
end

set(gca,'xtick',ticPos)
set(gca,'xticklabel',{'lh','rh'})
set(gca,'ytick',ticPos)
set(gca,'yticklabel',{'lh','rh'})

set(gca,'fontsize', 12);
% xlabel('ROIs', 'FontSize', 12)
% ylabel('ROIs', 'FontSize', 12)

% ti = get(gca,'TightInset');
% set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 2])

print('-dsvg','-r864',fullfile(plotDir,[fname '.svg']))
print('-dpng','-r864',fullfile(plotDir,[fname '.png']))
export_fig(fullfile(plotDir,[fname '.pdf']), '-transparent', '-r864');



end

