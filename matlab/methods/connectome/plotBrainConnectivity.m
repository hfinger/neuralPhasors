function [ figh ] = plotBrainConnectivity( nodeVals, edgeVals, nodeMinSize, nodeScale, nodeRange, edgeColor, edgeMin, edgeAlphas, surfaceVisibility, rescaleEdgeVals, resort )
%PLOTBRAINCONNECTIVITY Function that plots 66 Nodes in brain and connections between those nodes 
%   Input Parameters:                   
%       nodeVals - Vector of length 66 which is used to color the nodes
%       edgeVals - 66 x 66 matrix which is used to determine the
%                 linethickness of each edge
%       edgeColor - color specifier (string or triplet) used for edges
%       nodeMinSize - scalar that indicates minimum size of nodes
%       nodeScale - scalar, that determines how strongly the node size will
%                   be scaled according to nodeVals entries
%       nodeRange - min and max values for node coloring
%       edgeMin - Minimum value of edgeVals which is still shown in graph
%       edgeAlphas - vector of length size(edgeVals,3) that determines
%                    transparency of the edges on that dimension
%       surfaceVisibility - scalar between 0 and 1 that indicates 
%                           transparency of brain surface
%       rescaleEdgeVals - if true, edge values will be rescaled to give
%                         alternative edge thicknesses
%       resort - If true, rearange the order of nodeVals and edgeVals
%               (reverses order achieved by resortIdsMarlene)
%
%   Output:
%       figh - figure handle for graph

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
    if length(size(edgeVals)) < 3
        edgeVals = reshape(edgeVals, size(edgeVals,1), size(edgeVals,2), 1);
    end
    edgeVals = edgeVals(idx,:,:);
    edgeVals = edgeVals(:,idx,:);
end

%% plot nodes and edges on brain surface
 
metricVals = nodeVals;

% plot brain volume
figh = gcf;
clf;
set(figh,'Position', [200, 200, 700, 700])
set(gcf,'color','w');
clf;
hg=plot(g);
view(-9.8705,72.0327);
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
alpha(hg,surfaceVisibility);
hold on;

% colorvalues for nodes:
if isempty(nodeRange)
    cmin = min(metricVals);
    cmax = max(metricVals);
else
    cmin = nodeRange(1);
    cmax = nodeRange(2);
end
cmap=colormap('jet');
colormap(cmap);
yy=linspace(cmin,cmax,size(cmap,1));  % Generate range of color indices that map to cmap
cm = spline(yy,cmap',metricVals);     % Find interpolated colorvalue
cm(cm>1)=1;                           % Sometimes iterpolation gives values that are out of [0,1] range...
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
    conn = connRange ./ (maxConn - conn);
end

% set edge color and alpha
if isempty(edgeColor)
    edgeColor = hsv(size(conn,3));
elseif ischar(edgeColor)
    edgeColor = repmat(bitget(find('krgybmcw'==edgeColor)-1,1:3),size(conn,3),1);
elseif size(edgeColor,1) == 1
    edgeColor = repmat(edgeColor,size(conn,3),1);
end
if isempty(edgeAlphas)
    edgeAlphas = ones(size(conn,3));
end
edgeColor = horzcat(edgeColor,edgeAlphas');

% plot nodes and edges
for k=1:size(conn,1)  
    plot3(roi_center(k,1),roi_center(k,2),roi_center(k,3),'o','MarkerSize',nodeMinSize + abs(nodeVals(k)) * nodeScale,'color',cm(:,k), 'MarkerFaceColor', cm(:,k),'LineWidth', 2)
    for m=1:size(conn,2)
        for l=1:size(conn,3)
            if conn(m,k,l) > -inf && conn(m,k,l) ~= 0
                h=line([roi_center(k,1) roi_center(m,1)]',[roi_center(k,2) roi_center(m,2)]',[roi_center(k,3) roi_center(m,3)]','color',edgeColor(l,:),'LineWidth', conn(m,k,l));
            end
        end
    end  
end
hold off;

% general figure settings
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
OptionZ.FrameRate=10;
OptionZ.Duration=20;
OptionZ.Periodic=true;
view(0,40)


end

