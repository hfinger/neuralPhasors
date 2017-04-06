savefolder = '/net/store/nbp/projects/phasesim/workdir/20161026_fs_rois_subj3';
mkdir(savefolder);
subject = 3;
subjectStr = num2str(subject,'%02u');

data_path = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e';

g = gifti([data_path '/ca' subjectStr '_1_structcortex_8196.surf.gii']);

v1 = g.vertices(g.faces(:,1),:);
v2 = g.vertices(g.faces(:,2),:);
v3 = g.vertices(g.faces(:,3),:);

v1v2 = v2 - v1;
v2v3 = v3 - v2;

faceNorm = cross(v1v2, v2v3);
faceNorm = bsxfun(@rdivide,faceNorm,sqrt(sum(faceNorm.^2,2)));

faceCenter = (v1+v2+v3)/3;

% calc vertex normal vector by averaging over faces
vertexNorm = zeros(size(g.vertices));
for vid=1:size(g.vertices,1)
  fids = find(sum(g.faces==vid,2));
  vertexNorm(vid,:) = mean(faceNorm(fids,:),1);
end
vertexNorm = bsxfun(@rdivide,vertexNorm,sqrt(sum(vertexNorm.^2,2)));

%
fs_rois = load([data_path '/fs_rois/ca' subjectStr '_fs_rois.mat']);
regionmapping = importdata([data_path '/' subjectStr '_regionmapping.txt']);
lf = load([data_path '/ca' subjectStr '_lf_cortex.mat']);
SC=load(fullfile(data_path, '../..', 'Stuct_connectivity_for_Holger.mat'));

% roi_center = SC.roi_positions(cell2mat(SC.struct_labels_correspnding_pos),:);
roi_center = fs_rois.fs_rois;
% roi_center = lf.lf.pos;

%
numRois = 66;
roiNorm = zeros(numRois,3);
for k=1:numRois
  ids=find(regionmapping==k);
  roiNorm(k,:) = mean(vertexNorm(ids,:),1);
end
roiNormUnit = bsxfun(@rdivide,roiNorm,sqrt(sum(roiNorm.^2,2)));

% calculate roiIds of each Face
roiOfFace = regionmapping(g.faces);
roiOfFace2 = zeros(size(roiOfFace,1),1);

ids = roiOfFace(:,1)==roiOfFace(:,2);
roiOfFace2(ids) = roiOfFace(ids,1);
ids = roiOfFace(:,2)==roiOfFace(:,3);
roiOfFace2(ids) = roiOfFace(ids,2);
ids = roiOfFace(:,3)==roiOfFace(:,1);
roiOfFace2(ids) = roiOfFace(ids,3);

roiOfFace1 = zeros(size(roiOfFace,1),1);
ids = (roiOfFace(:,1)==roiOfFace(:,2)) & (roiOfFace(:,1)==roiOfFace(:,3));
roiOfFace1(ids) = roiOfFace(ids,1);

% color data for plots:
cdataPerROI = distinguishable_colors(33); %jet(67);
cdataPerROI=cdataPerROI(randperm(33),:);
cdataPerROI = [cdataPerROI; cdataPerROI; 0 0 0];

regionmappingPlusOne = regionmapping;
regionmappingPlusOne(regionmappingPlusOne==0) = 67;
cdataPerVertex.cdata = cdataPerROI(regionmappingPlusOne,:);


% load projection mat:
load(fullfile(data_path, '..', 'projectionMat/caAll_projectionMat.mat'))

electrodeLoc = importdata(fullfile(data_path, '..',['ca_electrodeLocations/ca' subjectStr '_EEGLocations.txt']));
electrodeLocLabels = electrodeLoc.textdata;
electrodeLoc = electrodeLoc.data;

fs_rois = load([data_path '/fs_rois/ca' subjectStr '_fs_rois.mat']);

reorderIds = cellfun(@(x) find(strcmp(electrodeLocLabels,x)), electrodeLabels);

dataSC=load(fullfile(data_path, '../..', 'dist_and_CI_controls_preprocessed.mat'));
avgSurfaceUniformity = load(fullfile(data_path,'../projectionMat/avgSurfaceUniformity.mat'));
SCmetrics = load(fullfile(data_path,'../../dist_and_CI_controls_network_metrics.mat'));
SCmetrics.perROI.surfUniformity = avgSurfaceUniformity.avgSurfaceUniformity;


views.posterior = [0 0];
views.leftLateral = [-90 0];
views.anterior = [180 0];
views.rightLateral = [90 0];
views.dorsal = [0 90];
views.ventral = [180 -90];

viewnames = fieldnames(views);

videoviews = [...
  0,40;...
  90,40;...
  180,40;...
  270,40;...
  360,-40;...
  450,-40;...
  540,-40;...
  630,-40;...
  720,-40;...
  ];





roiId=15;

load(fullfile(data_path, 'customCmap.mat'))

figh = sfigure(10);
clf;
set(figh,'Position', [400, 400, 500, 500])
set(gcf,'color','w');
% plot3(fs_rois.fs_rois(roiId,1),fs_rois.fs_rois(roiId,2),fs_rois.fs_rois(roiId,3),'ro');
% hold on;
% scatter3(electrodeLoc(reorderIds,1),electrodeLoc(reorderIds,2),electrodeLoc(reorderIds,3),10,projMatRoiPerVertex{subject}(:,roiId))

hg=plot(g);
alpha(hg,0.05);

hold on;
faceIds = find(sum(roiOfFace==roiId,2));
gROI=g;
gROI.faces = gROI.faces(faceIds,:);
%   hg=plot(gca,gROI,cdataPerVertex);
% hg=plot(gca,gROI);
% alpha(hg,0.25);
hold on;
dirFace=6;
dirCenter=20;

view(-9.8705,72.0327);

vertexIds = find(regionmapping==roiId);

plot3([g.vertices(vertexIds,1)]',[g.vertices(vertexIds,2)]',[g.vertices(vertexIds,3)]','o','Color','r')
plot3([g.vertices(vertexIds,1)+dirFace*vertexNorm(vertexIds,1) g.vertices(vertexIds,1)]',[g.vertices(vertexIds,2)+dirFace*vertexNorm(vertexIds,2) g.vertices(vertexIds,2)]',[g.vertices(vertexIds,3)+dirFace*vertexNorm(vertexIds,3) g.vertices(vertexIds,3)]','Color','r')
plot3(roi_center(roiId,1),roi_center(roiId,2),roi_center(roiId,3),'o','Color','k','LineWidth',3)
plot3([roi_center(roiId,1)+dirCenter*roiNormUnit(roiId,1) roi_center(roiId,1)]',[roi_center(roiId,2)+dirCenter*roiNormUnit(roiId,2) roi_center(roiId,2)]',[roi_center(roiId,3)+dirCenter*roiNormUnit(roiId,3) roi_center(roiId,3)]','Color','k','LineWidth',3)

%   scatter3(electrodeLoc(reorderIds,1),electrodeLoc(reorderIds,2),electrodeLoc(reorderIds,3),50,projMatRoiPerVertex{subject}(:,roiId))
% scatter3(electrodeLoc(reorderIds,1),electrodeLoc(reorderIds,2),electrodeLoc(reorderIds,3),50,projMatRoi{subject}(:,roiId))

colormap(cmap)

c=projMatRoi{subject}(:,roiId);
cmin=min(c);
cmax=max(c);
cmap=colormap(cmap);                      % Set colormap
yy=linspace(cmin,cmax,size(cmap,1));  % Generate range of color indices that map to cmap
cm = spline(yy,cmap',c);                  % Find interpolated colorvalue
cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;
caxis([cmin cmax])

for k=1:63
  plot3(electrodeLoc(reorderIds(k),1),electrodeLoc(reorderIds(k),2),electrodeLoc(reorderIds(k),3),'o','MarkerSize',10,'color',[cm(:,k)],'LineWidth', 3)
end


colorbar
axis off;

set(gca,'FontSize',12)
% handle=title({['ROI ' num2str(roiId)]; SC.struct_labels{roiId}});
% set(handle,'Position',[0,60,0]);

% export_fig(fullfile(savefolder,['ROI_lf_' num2str(roiId,'%02u') '.pdf']), '-opengl');

% export figure as png without colorbar
colorbar('off')
export_fig(fullfile(savefolder,['ROI_lf_' num2str(roiId,'%02u') '.png']), '-r200');

% export only colorbar as pdf
figh = sfigure(11);
clf;
set(gcf,'color','w');
set(figh,'Position', [100, 100, 200, 200])
colormap(cmap)
caxis([cmin cmax])
set(gca,'FontSize',12)
c = colorbar;
x1=get(gca,'position');
x=get(c,'Position');
x(3)=0.03;
set(c,'Position',x)
x1(1)=-5+x1(1);
set(gca,'position',x1)
set(gca,'fontsize', 12);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 2])
export_fig(fullfile(savefolder,['ROI_lf_' num2str(roiId,'%02u') '_colorbar.pdf']));
