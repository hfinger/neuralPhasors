function cortexTopology(SC, SCMetrics,aaa)
%% load surface and ROI 3D data:
datapath = '/net/store/nbp/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e';

subject = 3;
subjectStr = num2str(subject,'%02u');

g = gifti(fullfile(datapath,['ca' subjectStr '_1_structcortex_8196.surf.gii']));
fs_rois = load(fullfile(datapath,['fs_rois/ca' subjectStr '_fs_rois.mat']));
%regionmapping = importdata(fullfile(datapath,[subjectStr '_regionmapping.txt']));
%lf = load(fullfile(datapath,['ca' subjectStr '_lf_cortex.mat']));

roi_center = fs_rois.fs_rois;

%paths = dataPaths();
%dist = load(fullfile(paths.databases,'SC_Bastian','icoh_all_lcmvhilbertrest_20140807.mat'),'pos');
%dist = squeeze(sqrt(sum(bsxfun(@minus,dist.pos,permute(dist.pos,[3 2 1])).^2,2)));
%avgSurfaceUniformity = load([paths.databases '/SC_Bastian/surfaces/projectionMat/avgSurfaceUniformity.mat']);

%% plot figure

figh = sfigure();
set(figh,'Position', [200, 200, 700, 700])
set(gcf,'color','w');
clf;
hg=plot(g);
view(-9.8705,72.0327);
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
alpha(hg,0.15);
hold on;

% color values:
cmin = min(SCMetrics.perROI.clustCoef);
cmax = max(SCMetrics.perROI.clustCoef);
cmap=colormap('jet');
colormap(cmap);
yy=linspace(cmin,cmax,size(cmap,1));                                        % Generate range of color indices that map to cmap
cm = spline(yy,cmap',SCMetrics.perROI.clustCoef);                           % Find interpolated colorvalue 
cm(cm>1)=1;                                                                 % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;     
caxis([cmin cmax])

% plot only the 600 max values
SC=log(SC);
SC(isnan(SC))=-Inf;
connSorted = sort(SC(:));
connmin = connSorted(end-600);
SC(SC<connmin) = -Inf;
connmax=max(SC(:));

aaa = (aaa==2);

for k=1:66
    
    if(aaa(k))>0
        plot3(roi_center(k,1),roi_center(k,2),roi_center(k,3),'o','MarkerSize',13,'color',cm(:,k), 'MarkerFaceColor', cm(:,k),'LineWidth', 2)
        for m=1:k-1
            if 5*(SC(m,k)-connmin)/(connmax-connmin)>0
                h=line([roi_center(k,1) roi_center(m,1)]',[roi_center(k,2) roi_center(m,2)]',[roi_center(k,3) roi_center(m,3)]',...
                    'color','k','LineWidth',5*(SC(m,k)-connmin)/(connmax-connmin));
            end
        end
    end
end

hold off;
colorbar
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])