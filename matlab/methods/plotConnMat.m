function plotConnMat( connMat, plotDir, fname, clims )
%PLOTCONNMAT Summary of this function goes here
%   Detailed explanation goes here

figh = sfigure(1);
clf;
set(figh,'Position', [100, 100, 200, 200])
imagesc(connMat)
axis image;

% title(['average subject SC'])
% colormap hot;

c = colorbar;
if nargin>3
  set(gca,'clim',clims);
end
x1=get(gca,'position');
x=get(c,'Position');
x(3)=0.03;
set(c,'Position',x)
x1(1)=0.08+x1(1);
set(gca,'position',x1)

set(gca,'xtick',[16 49])
set(gca,'xticklabel',{'lh','rh'})
set(gca,'ytick',[16 49])
set(gca,'yticklabel',{'lh','rh'})

set(gca,'fontsize', 12);
xlabel('ROIs', 'FontSize', 12)
ylabel('ROIs', 'FontSize', 12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 2])

print('-dpng','-r864',fullfile(plotDir,[fname '.png']))
export_fig(fullfile(plotDir,[fname '.pdf']), '-transparent', '-r864');



end

