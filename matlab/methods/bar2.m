function bar2( data, legend1, legend2, plotDir, fname, xlimval, ylimval )
%PLOTBAR2 Summary of this function goes here
%   Detailed explanation goes here

figh = sfigure(1);
clf;
set(figh,'Position', [100, 100, 450, 200])
bar(data)

set(gca,'fontsize', 12);
legend(legend2, 'FontSize', 12, 'Location','NorthEastOutside')
set(gca,'XTickLabel',legend1, 'FontSize', 12)
ylabel('correlation', 'FontSize', 12)

if nargin>5 && ~isempty(xlimval)
  xlim(xlimval)
end
if nargin>6 && ~isempty(ylimval)
  ylim(ylimval)
end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 2])

print('-dpng','-r72',fullfile(plotDir,[fname '.png']))
export_fig(fullfile(plotDir,[fname '.pdf']), '-transparent', '-r72');

end

