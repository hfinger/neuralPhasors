function cortexTopology(SCMetrics,SC,sub)

% SC,SCMetrics: /net/store/nbp/projects/phasesim/workdir/pebel/20150414_SAR_Metrics/ConnectomeMetrics
%               /net/store/nbp/projects/phasesim/workdir/pebel/20150414_SAR_Metrics/temp_ConnectomeMetrics
% with param k: /net/store/nbp/projects/phasesim/workdir/pebel/20150414_SAR_Metrics/ConnectomeSim

%% get data and call function
paths = dataPaths();
datapath = strcat(paths.databases,'/SC_Bastian/surfaces/wetransfer-b16a3e');
plotsTo =  strcat(paths.localTempDir,'/Results/','Topography/');
subject = 3; subjectStr = num2str(subject,'%02u');

load(strcat(paths.databases ,'/SC_Bastian/SC_labelsPairs.mat'))
%load('C:\Users\PWJEbel\Desktop\USB goes HERE\this\databases\SC_Bastian\SC_labelsPairs.mat')

s.g = gifti(fullfile(datapath,['ca' subjectStr '_1_structcortex_8196.surf.gii']));
s.fs_rois = load(fullfile(datapath,['fs_rois/ca' subjectStr '_fs_rois.mat']));
s.roi_center = s.fs_rois.fs_rois;

if nargin < 2
  load(fullfile(paths.databases,'SC_Bastian','dti_20141209_preprocessed.mat'));
  SC = avg_ci;
  SC(isnan(SC)) = 0;
  SC = SC + SC';
  SC = normGraph(SC, avg_roi_size, 'ROIprd', true, 0);
  % sub can be used to specify a subset of nodes to be plotted (communalities, rich club, outlier metrics ...)
end

if nargin < 3
  sub = ones(size(SC,1));
end

%paths = dataPaths();
%dist = load(fullfile(paths.databases,'SC_Bastian','icoh_all_lcmvhilbertrest_20140807.mat'),'pos');
%dist = squeeze(sqrt(sum(bsxfun(@minus,dist.pos,permute(dist.pos,[3 2 1])).^2,2)));
%avgSurfaceUniformity = load([paths.databases '/SC_Bastian/surfaces/projectionMat/avgSurfaceUniformity.mat']);

for i= [2,5,6,7,8,18] 
  plotter(plotsTo, SC,SCMetrics, i,s, sub, false, SC_labelsPairs)
  plotter(plotsTo, SC,SCMetrics, i,s, sub, true, SC_labelsPairs)
end

end

%% function that plots cortical surface, nodes and edges
function plotter(path, SC,SCMetrics, i, s, sub, asym, SC_labelsPairs)

if(asym) path = strcat(path,'Asymmetry/'); end
if ~exist(path,'dir')
  mkdir(path);
end

ut = triu(true(size(SC)),+1)+triu(true(size(SC)),+1)';                      % de-select main diagonal
inter = logical(zeros(size(SC)));                                           % mask for interhemispheric connections
inter(1:length(SC)/2,length(SC)/2+1:end) = true;
interRL = logical(inter);
interLR = logical(inter');
inter = logical(interRL + interLR);
intra = logical(logical(ut - inter));                                       % mask for intrahemispheric connections
intraLL = and(intra, vertcat(ones(length(SC)/2,length(SC)),zeros(length(SC)/2,length(SC))));
intraRR = and(intra, vertcat(zeros(length(SC)/2,length(SC)),ones(length(SC)/2,length(SC))));

if i == 1 % little hack bc data field was inserted later on
  %label = 'shortestPaths';
  label = SCMetrics.lables.perROI{i};
else
  label = SCMetrics.lables.perROI{i};
end

metr = SCMetrics.perROI.(label);                                            % get metrics to be color-coded, select from SCMetrics.perROI

if(asym)                                                                    % calculate pairwise differences in metrics: LH - RH
  diff = metr(1:length(SC)/2) - metr(length(SC)/2 + 1:length(SC));          % values > 0: left-dominance, values < 0: right-dominance
  norm = mean([metr(1:length(SC)/2), metr(length(SC)/2 + 1:length(SC))],2);
  metr = 100* diff ./ norm;
  metr(isnan(metr),:) = 0;
  % remove outliers? 200% deviation bc one value is zero (other might be close to zero)
  metr(find(metr==-200,numel(metr)),:) = 0;
  metr(find(metr==+200,numel(metr)),:) = 0;
  % paired differences outlier detection:
  figure();subplot(1,2,1);boxplot(metr);
  subplot(1,2,2); qqplot(metr); title(label);                         
  print(strcat(path,strcat('/box_',label)),'-dpng');
  % Grubb's test (req. normality), Dixon's Q test (robust to normality assumption in small samples)
  
  % spatial outlier detection: Kou, Lu, Chen (2006)
  metr = vertcat(metr,metr);
end

figh = figure();
set(figh,'Position', [200, 200, 700, 700])
set(gcf,'color','w');
clf;
hg=plot(s.g);
view(-9.8705,72.0327);
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
alpha(hg,0.15);
hold on;

% color values:
cmin = min(metr);
cmax = max(metr);
cmap=colormap(jet);
if(asym) cmap=colormap(b2r(cmin,cmax)); end
colormap(cmap);
yy=linspace(cmin,cmax,size(cmap,1));                                        % Generate range of color indices that map to cmap
cm = spline(yy,cmap',metr);                                                 % Find interpolated colorvalue
cm(cm>1)=1;                                                                 % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;
caxis([cmin cmax])

% plot only the 400 max values
SC=log(SC);
SC(isnan(SC))=-Inf;
connSorted = sort(SC(:));
connmin = connSorted(end-400);
SC(SC<connmin) = -Inf;
connmax=max(SC(:));

for k=1:size(SC,1)
  if(sub(k))>0  % if ROI k is in selected subset, do:
    plot3(s.roi_center(k,1),s.roi_center(k,2),s.roi_center(k,3),'o','MarkerSize',...
          13,'color',cm(:,k), 'MarkerFaceColor', cm(:,k),'LineWidth', 2)
    for m=1:k-1
      if 5*(SC(m,k)-connmin)/(connmax-connmin)>0
        h=line([s.roi_center(k,1) s.roi_center(m,1)]',[s.roi_center(k,2) s.roi_center(m,2)]',...
          [s.roi_center(k,3) s.roi_center(m,3)]','color','k','LineWidth',5*(SC(m,k)-connmin)/(connmax-connmin));
      end
    end
  end
end

hold off;
handle = colorbar;
if(asym) ylabel(handle, '% of directed deviation from pairwise mean'); end;          % perc. asymm. = (LH-RH)./mean(LH,RH)
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
str = label;
if(asym) str = strcat(str, ' (asymmetry)'); end
text(50, 50, 0, str)
set(gcf,'PaperPositionMode','auto')
print(strcat(path,label),'-dpng', '-r0');


if(asym)
    metr = SCMetrics.perROI.(label);
    diff = metr(1:length(SC)/2) - metr(length(SC)/2 + 1:length(SC));          % values > 0: left-dominance, values < 0: right-dominance
    norm = mean([metr(1:length(SC)/2), metr(length(SC)/2 + 1:length(SC))],2);
    metr = 100* diff ./ norm;
    metr(isnan(metr),:) = 0;
    % remove outliers? 200% deviation bc one value is zero (other might be close to zero)
    metr(find(metr==-200,numel(metr)),:) = 0;
    metr(find(metr==+200,numel(metr)),:) = 0;
    
    % get bar plots
    [metrSort, I] = sort(abs(metr), 'descend'); % take totals, color-code red(>0, LH) and blue(<0, RH), gray for < 10% dev
    metrSortSign = metr(I);
    metrColor = 0.5*ones(length(metrSort),3);
    %metrColor(:,1) = 1;
    metrColor(metrSortSign>0,:) = repmat([1 0 0], sum(metrSortSign>0,1),1);
    metrColor(metrSortSign<0,:) = repmat([0 0 1], sum(metrSortSign<0,1),1);
    
    fHand = figure;
    aHand = axes('parent', fHand);
    hold(aHand, 'on')
    for j = 1:numel(metrSort)
        bar(j, metrSort(j), 'parent', aHand, 'facecolor', metrColor(j,:));
    end
    
    title(label)
    ylabel('abs. percent of deviation')
    set(gca,'XTick',1:length(metrSort))
    set(gca,'XTickLabel',SC_labelsPairs(I))%{'date','date','date','date'})
    rotateXLabels(aHand, 35)%, varargin )
    set(fHand, 'Position', [0 0 900 300])
    %print(strcat(pwd,SCMetrics.lables.perROI{i},'Bar'),'-dpng', '-r0')
    img = getframe(gcf);
    imwrite(img.cdata, [strcat(path,label,'Bar'), '.png']);
end

end

