function [ hg ] = plotBrainConnectivity( nodeVals, edgeVals, edgeMin, resort, rescaleEdgeVals )
%PLOTBRAINCONNECTIVITY Function that plots 66 Nodes in brain and connections between those nodes 
%   Input Parameters:                   
%       nodeVals - Vector of length 66 which is used to color the nodes
%       edgeVals - 66 x 66 matrix which is used to determine the
%                 linethickness of each edge
%       edgeMin - Minimum value of edgeVals which is still shown in graph
%       resort - If true, rearange the order of nodeVals and edgeVals
%               (reverses order achieved by resortIdsMarlene)
%
%   Output:
%       hg - figure handle for graph

%% load brain surface
subject = 3;
savefolder = strcat('Coherence_Stimulated_subj',num2str(subject));
mkdir(savefolder);
subjectStr = num2str(subject,'%02u');

addpath('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e')
g = gifti(['ca' subjectStr '_1_structcortex_8196.surf.gii']);
fs_rois = load(['fs_rois/ca' subjectStr '_fs_rois.mat']);
roi_center = fs_rois.fs_rois;

%% resort data if necessary
if resort
    path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
    load(path_ResortIDs)
    resortIDs = [resortIdsMarlene, resortIdsMarlene + 33];
    idx = 1:length(resortIDs);
    for i=idx
        idx(resortIDs(i)) = i;
    end
    nodeVals = nodeVals(idx);
    edgeVals = edgeVals(idx,:);
    edgeVals = edgeVals(:,idx);
end

%% plot nodes and edges on brain surface
 
metricVals = nodeVals;

% plot brain volume
figh = sfigure();
clf;
set(figh,'Position', [200, 200, 700, 700])
set(gcf,'color','w');
clf;
hg=plot(g);
view(-9.8705,72.0327);
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
alpha(hg,0.15);
hold on;

% colorvalues for nodes:
cmin = min(metricVals);
cmax = max(metricVals);
cmap=colormap('jet');
colormap(cmap);
yy=linspace(cmin,cmax,size(cmap,1));  % Generate range of color indices that map to cmap
cm = spline(yy,cmap',metricVals);                  % Find interpolated colorvalue
cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;
caxis([cmin cmax])

% scale edge values
conn=edgeVals;
conn(conn<edgeMin) = 0;
if rescaleEdgeVals
    conn = log(conn);
    minConn = min(conn(conn > -inf));
    maxConn = 0;
    connRange = maxConn - minConn;
    newConn = connRange ./ (maxConn - conn);
else
    newConn = conn;
end
% plot nodes and edges
for k=1:66  
    plot3(roi_center(k,1),roi_center(k,2),roi_center(k,3),'o','MarkerSize',13,'color',cm(:,k), 'MarkerFaceColor', cm(:,k),'LineWidth', 2)
    for m=1:66
        if conn(m,k) > -inf && conn(m,k) ~= 0
            h=line([roi_center(k,1) roi_center(m,1)]',[roi_center(k,2) roi_center(m,2)]',[roi_center(k,3) roi_center(m,3)]','color','r','LineWidth', newConn(m,k));
        end
    end  
end
hold off;
%   colorbar
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])

OptionZ.FrameRate=10;
OptionZ.Duration=20;
OptionZ.Periodic=true;
view(0,40)

% export only colorbar as pdf
%figh = sfigure();
%clf;
%set(gcf,'color','w');
%set(figh,'Position', [100, 100, 200, 200])
%colormap(cmap)
%caxis([cmin cmax])
%set(gca,'FontSize',12)
%c = colorbar;
%x1=get(gca,'position');
%x=get(c,'Position');
%x(3)=0.03;
%set(c,'Position',x)
%x1(1)=-5+x1(1);
%set(gca,'position',x1)
%set(gca,'fontsize', 12);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 2])
%export_fig(fullfile(savefolder,['SC_' fnames{fid} '_colorbar.pdf']));

end

