function [ figh ] = plotBrainConnectivityNew( nodeVals, edgeVals, nodeMinSize, nodeScale, nodeRange, edgeColor, edgeMin, edgeAlphas, surfaceVisibility, rescaleEdgeVals, resort, cmap )
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

dpaths = dataPaths();
matlab_path = dpaths.sge_pathToAddScriptPaths;

% load brain surface and plot it:
g = gifti(fullfile(matlab_path, 'methods/connectome/brain_plot_data/surfaceFsData/ca01_1_structcortex_8196.surf.gii'));

% find all vertices on right hemisphere:
vertices_left = find(g.vertices(:,1)<0);
keep_faces = ismember(g.faces(:,1), vertices_left) & ismember(g.faces(:,2), vertices_left) & ismember(g.faces(:,3), vertices_left);
g.faces = g.faces(keep_faces,:);


%% load tractography data:
freesurfer_roi_ids = load(fullfile(matlab_path, 'methods/connectome/brain_plot_data/tractographyData/freesurfer_roi_ids.mat'));
freesurfer_roi_ids = freesurfer_roi_ids.freesurfer_roi_ids;

voxel_coords = load(fullfile(matlab_path, 'methods/connectome/brain_plot_data/tractographyData/voxel_coords.mat'));
voxel_coords = voxel_coords.voxel_coords;

%% calculate center of each roi:
roi_center_coords = zeros(66,3);
for k=1:66
    roi_center_coords(k,:) = mean(voxel_coords(freesurfer_roi_ids==k,:),1);
end

%% load nifti header to calculate transform matrix (from Tractography Space to Freesurfer Space):
header = cbiReadNiftiHeader(fullfile(matlab_path, 'methods/connectome/brain_plot_data/surfaceFsData/compl_fs_mask.nii'));
srow_mat = [header.srow_x; header.srow_y; header.srow_z];

%% load voxel coordinates in tractography space:
roi_center = [roi_center_coords ones(size(roi_center_coords,1),1)] * srow_mat';

%% resort data if necessary
if resort
    path_ResortIDs = fullfile(matlab_path, 'methods/connectome/brain_plot_data/sources_plot_order.mat');
    resortIdsMarlene = load(path_ResortIDs);
    resortIDs = resortIdsMarlene.sort_ids;
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
    
    if length(size(edgeColor)) < 3
        edgeColor = reshape(edgeColor, size(edgeColor,1), size(edgeColor,2), 1);
    end
    edgeColor = edgeColor(idx,:,:);
    edgeColor = edgeColor(:,idx,:);
end

%% plot nodes and edges on brain surface
 
metricVals = nodeVals;

% plot brain volume
figh = gcf;
clf;
set(figh,'Position', [200, 200, 1200, 900])
set(gcf,'color','w');
clf;
hg=plot(g);
alpha(hg,0.1);
set(gca,'xtick',[],'ytick',[])
view(-9.8705,72.0327);
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
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
colormap(cmap);
yy=linspace(cmin,cmax,size(cmap,1));  % Generate range of color indices that map to cmap
cm = spline(yy,cmap',metricVals);     % Find interpolated colorvalue
cm(cm>1)=1;                           % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;
caxis([cmin cmax])

% scale edge values
conn=edgeVals;
conn(conn<edgeMin) = 0;
if rescaleEdgeVals==2
    minConn = min(conn(conn > -inf));
    maxConn = max(conn(conn < inf));
    connRange = maxConn - minConn;
    conn = 5*(conn - minConn) / connRange;
elseif rescaleEdgeVals
    conn = log(conn);
    minConn = min(conn(conn > -inf));
    maxConn = 0;
    connRange = maxConn - minConn;
    conn = connRange ./ (maxConn - conn);
end

% colorvalues for edges:
cminedge = min(edgeColor(:));
cmaxedge = max(edgeColor(:));
colormap(cmap);
yy=linspace(cminedge,cmaxedge,size(cmap,1));  % Generate range of color indices that map to cmap
cmedge = spline(yy,cmap',edgeColor);     % Find interpolated colorvalue
cmedge(cmedge>1)=1;                           % Sometimes iterpolation gives values that are out of [0,1] range...
cmedge(cmedge<0)=0;

% plot nodes and edges
for k=1:size(conn,1)
    current_marker_size = nodeMinSize + abs(nodeVals(k)) * nodeScale;
    if current_marker_size > 0
        plot3(roi_center(k,1),roi_center(k,2),roi_center(k,3),'o','MarkerSize',current_marker_size,'color',cm(:,k), 'MarkerFaceColor', cm(:,k),'LineWidth', 2)
    end
    for m=1:size(conn,2)
        if conn(m,k) > -inf && conn(m,k) ~= 0
            h=line([roi_center(k,1) roi_center(m,1)]',[roi_center(k,2) roi_center(m,2)]',[roi_center(k,3) roi_center(m,3)]','color',cmedge(:,m,k),'LineWidth', conn(m,k));
        end
    end
end

hold off;

% general figure settings
set(gca,'FontSize',12)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8])
OptionZ.FrameRate=10;
OptionZ.Duration=20;
OptionZ.Periodic=true;
view(0,40)


end

