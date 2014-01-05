function plot4D(data, dimName, dimLabels, caption, filename, plotDir, permutePlotDims)

if nargin<7
  permutePlotDims = [];
end

%% plot:
siz = size(data);

if ~isempty(permutePlotDims)
  data = permute(data,permutePlotDims);
  while length(dimName)<length(permutePlotDims)
    dimName{end+1} = '';
    dimLabels{end+1} = {};
  end
  dimName = dimName(permutePlotDims);
  dimLabels = dimLabels(permutePlotDims);
end

if length(siz) > 2
  if length(siz) == 3
    %% shift third dimension to forth for horizonal plot
    data = permute(data,[1 2 4 3]);
    dimName{4} = '';
    dimLabels{4} = {};
    dimName = dimName([1 2 4 3]);
    dimLabels = dimLabels([1 2 4 3]);
  end
  
  %           if length(siz) == 3
  %
  %             m=round(sqrt(siz(3)));
  %             n=ceil(siz(3)/m);
  %
  %             sfigure(1); clf;
  %             for k=1:siz(3)
  %               subplot(m,n,k);
  %               imagesc(data(:,:,k));
  %               title(caption)
  %               colormap hot;
  %               colorbar;
  %               if length(dimName)>1
  %                 xlabel(dimName{2})
  %                 set(gca,'XTick',1:size(data,2));
  %                 set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),dimLabels{2},'UniformOutput',false))
  %               end
  %               ylabel(dimName{1})
  %               set(gca,'YTick',1:size(data,1));
  %               set(gca,'YTickLabel',dimLabels{1})
  %             end
  %             set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
  %             print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
  %
  %           elseif length(siz) == 4
  
  minc = min(data(:));
  maxc = max(data(:));
  
  h=sfigure(1);
  clf
  p = panel();
  p.pack(size(data,3), size(data,4));
  p.de.margin = 2;
  p.margin = [25 25 20 20];
  p.fontsize = 8;
  p.title(caption);
  for m = 1:size(data,3)
    for n = 1:size(data,4)
      p(m, n).select();
      imagesc(data(:,:,m,n));
      
      %                 [nr,nc] = size(data(:,:,m,n));
      %                 pcolor([data(:,:,m,n) nan(nr,1); nan(1,nc+1)]);
      %                 shading flat;
      %                 set(gca, 'ydir', 'reverse');
      
      colormap autumn;
      if isnan(minc) || isnan(maxc) || minc==maxc || isinf(maxc) || isinf(minc)
        disp('cannot set caxis limits')
      else
        set(gca,'clim',[minc maxc]);
      end
      
      axis image;
      
      if n==size(data,4)
        originalSize = get(gca, 'Position');
        colorbar;
        set(gca, 'Position', originalSize);
      end
      
      if m==size(data,3)
        if n==round(size(data,4)/2)
          %% with outer xlabel
          p(m,n).xlabel({dimName{2},[],['{\bf' num2str(dimLabels{4}{n}) '}'],['{\bf' dimName{4} '}']});
        else
          %% only inner xlabel
          p(m,n).xlabel({dimName{2},[],['{\bf' num2str(dimLabels{4}{n}) '}']});
        end
        %set(gca,'XTick',1:size(data,1));
        xticks = get(gca,'XTick');
        xticks(xticks>length(dimLabels{2})) = [];
        xticks = xticks(mod(xticks,1)==0);
        set(gca,'XTick',xticks);
        if ischar(dimLabels{2}{1})
          set(gca,'XTickLabel',dimLabels{2}(xticks))
        else
          set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),dimLabels{2}(xticks),'UniformOutput',false))
        end
        %set(gca,'XTickMode','auto')
      else
        set(gca, 'xticklabel', {});
      end
      if n==1
        if isempty(dimLabels{3})
          p(m,n).ylabel(dimName{1});
        else
          if m==round(size(data,3)/2)
            %% with outer ylabel
            p(m,n).ylabel({['{\bf' dimName{3} '}'],['{\bf' num2str(dimLabels{3}{m}) '}'],[],dimName{1}});
          else
            %% only inner ylabel
            p(m,n).ylabel({['{\bf' num2str(dimLabels{3}{m}) '}'],[],dimName{1}});
          end
        end
        yticks = get(gca,'YTick');
        yticks(yticks>length(dimLabels{1})) = [];
        yticks = yticks(mod(yticks,1)==0);
        set(gca,'YTick',yticks);
        set(gca,'YTickLabel',dimLabels{1}(yticks))
      else
        set(gca, 'yticklabel', {});
      end
      
    end
  end
  
  
  %             set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
  %             print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
  %             export_fig(fullfile(plotDir,[filename '.pdf']))
  %             p.export(fullfile(plotDir,[filename '.pdf']), '-a1.4', '-rp');
  p.export(fullfile(plotDir,[filename '.pdf']), '-rp','-w300','-h100');
  saveas(h,fullfile(plotDir,[filename '.fig']));
  
  %           else
  %             return;
  %           end
elseif siz(2) == 1 && siz(1) == 1
  %% return because only scalar
  return;
elseif siz(2) == 1
  sfigure(1); clf;
  plot(data);
  title(caption)
  xlabel(dimName{1})
  ylabel(caption)
  set(gca,'XTick',1:size(data,1));
  set(gca,'XTickLabel',dimLabels{1})
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
  %           print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
  export_fig(fullfile(plotDir,[filename '.pdf']), '-transparent', '-r72');
else
  h=sfigure(1);clf;
  imagesc(data);
  title(caption)
  colormap hot;
  colorbar;
  if length(dimName)>1
    xlabel(dimName{2})
    set(gca,'XTick',1:size(data,2));
    set(gca,'XTickLabel',cellfun(@(x) num2str(x,'%.3g'),dimLabels{2},'UniformOutput',false))
  end
  ylabel(dimName{1})
  set(gca,'YTick',1:size(data,1));
  set(gca,'YTickLabel',dimLabels{1})
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
  %           print('-dpng','-r72',fullfile(plotDir,[filename '.png']))
  saveas(h,fullfile(plotDir,[filename '.fig']));
  export_fig(fullfile(plotDir,[filename '.pdf']), '-transparent', '-r72');
end

end
