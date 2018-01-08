
%% load surface
datapath = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e';
g = gifti(fullfile(datapath,['ca01_1_structcortex_8196.surf.gii']));

%% load nifti header to calculate transform matrix
header = cbiReadNiftiHeader('/net/store/nbp/projects/phasesim/databases/SC_Bastian/mosaic_2015_02_18/compl_fs_mask.nii');
srow_mat = [header.srow_x; header.srow_y; header.srow_z];

%% load gray matter coordinates
new_tract_space = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150423gridjob/new_tract_space.mat');
new_tract_space = new_tract_space.new_tract_space;

%% load selection ids to remove zero rows:
useVoxelIdx1 = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/useVoxelIdx1.mat');
useVoxelIdx1 = useVoxelIdx1.useVoxelIdx;
useVoxelIdx2 = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/useVoxelIdx2.mat');
useVoxelIdx2 = useVoxelIdx2.useVoxelIdx;
useVoxelIdx3 = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/useVoxelIdx3.mat');
useVoxelIdx3 = useVoxelIdx3.useVoxelIdx;
useVoxelIdx4 = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/useVoxelIdx4.mat');
useVoxelIdx4 = useVoxelIdx4.useVoxelIdx;
useVoxelIdx5 = load('/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/useVoxelIdx5.mat');
useVoxelIdx5 = useVoxelIdx5.useVoxelIdx;
useVoxelIdx = useVoxelIdx1(useVoxelIdx2(useVoxelIdx3(useVoxelIdx4(useVoxelIdx5))));

%% calculate for the used voxels the correct coordinates:
usedVoxelCoord = new_tract_space(useVoxelIdx,:);
usedVoxelCoord = [usedVoxelCoord ones(size(usedVoxelCoord,1),1)] * srow_mat';

%% select the closest gray matter voxel for each gifti vertex
distBetweenVoxAndVertex = zeros(size(g.vertices,1),1);
voxIdOfVertex = zeros(size(g.vertices,1),1);
for i=1:size(g.vertices,1)
  [distBetweenVoxAndVertex(i), voxIdOfVertex(i)] = min( sum( bsxfun(@minus, g.vertices(i,:), usedVoxelCoord).^2 ,2) );
end

%% select the closest gray matter voxel for each gifti face
faceCenterCoord = mean(cat( 3, g.vertices(g.faces(:,1),:), g.vertices(g.faces(:,2),:), g.vertices(g.faces(:,3),:) ),3);
distBetweenVoxAndFace = zeros(size(faceCenterCoord,1),1);
voxIdOfFace = zeros(size(faceCenterCoord,1),1);
for i=1:size(g.faces,1)
  [distBetweenVoxAndFace(i), voxIdOfFace(i)] = min( sum( bsxfun(@minus, faceCenterCoord(i,:), usedVoxelCoord).^2 ,2) );
end


%%
cdataPerROI = distinguishable_colors(120);
cdataPerROI = cdataPerROI(randperm(120),:);
az=0;
el=25;

%%
doCreateVideo = false;

%%
if doCreateVideo
  daObj=VideoWriter('topographyFsCos'); %for default video format. 
  daObj.FrameRate=10;
  open(daObj);
end

for i=68:100
  
  %%
  tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fullcos/clusterCenter' num2str(i-1) '.mat']);
  voxelIndexCoordsPrevious = tmp.clusterCenter;
    
  tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fullcos/clusterCenter' num2str(i) '.mat']);
  voxelIndexCoords = tmp.clusterCenter;
  tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fullcos/fullcosCluster' num2str(i) '.mat']);
  clusterIds = tmp.Cluster;
  tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fullcos/clusterConnmat' num2str(i) '.mat']);
  clusterConnmat = tmp.clusterConnmat;
  
  indexOfSplitCluster = find(sum(abs(voxelIndexCoordsPrevious-voxelIndexCoords(1:end-1,:)),2));
  
  clusterOfVertex = clusterIds(voxIdOfVertex);
  clusterOfFace = clusterIds(voxIdOfFace);
  
%   cdataPerVertex.cdata = cdataPerROI(clusterOfVertex,:);

  %%
  figh = figure(3);
  clf;
  
  cdataPerVertex.cdata = repmat([0.5 0.5 0.5],size(g.vertices,1),1);
%   hg=plot(g,cdataPerVertex);
%   alpha(hg,0.1);
%   shading flat
%   view(-9.8705,72.0327);
  
  
  
  
  twoColors = distinguishable_colors(2);
  faceIds1 = find(clusterOfFace==indexOfSplitCluster);
  faceIds2 = find(clusterOfFace==i);
  vertices1 = g.faces(faceIds1,:);
  vertices2 = g.faces(faceIds2,:);
  vertices1 = unique(vertices1(:));
  vertices2 = unique(vertices2(:));
  gROI=g;
%   gROI.faces = gROI.faces([faceIds1; faceIds2],:);
  cdataPerVertex.cdata = repmat([0.5 0.5 0.5],size(g.vertices,1),1);
  cdataPerVertex.cdata(vertices1,:) = repmat([0 0 1],length(vertices1),1);
  cdataPerVertex.cdata(vertices2,:) = repmat([1 0 0],length(vertices2),1);
  hg=plot(gROI,cdataPerVertex);
  alphaVal = 0.7*ones(size(g.vertices,1),1);
  alphaVal([vertices1; vertices2]) = 1;
  shading flat
  set(hg,'FaceAlpha','interp')
  set(hg,'FaceVertexAlphaData',alphaVal)
%   alpha(hg,0.5);
  set(gca,'xtick',[],'ytick',[])
  
  hold on;

  clusterCenterCoord = [voxelIndexCoords ones(size(voxelIndexCoords,1),1)] * srow_mat';

%   figh = sfigure(3);
%   clf;
%   set(figh,'Position', [200, 200, 700, 700])
%   set(gcf,'color','w');
%   clf;
%   hg=plot(g);
%   view(-9.8705,72.0327);
%   set(gca,'FontSize',12)
%   set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4])
%   alpha(hg,0.7);
  hold on;

  plot3(clusterCenterCoord(:,1),clusterCenterCoord(:,2),clusterCenterCoord(:,3),'o','MarkerFaceColor',[0 0 0],'markersize',6)

  
  
%   conn=log(conn);
%   conn(isnan(conn))=-Inf;
%   connSorted = sort(conn(:));
%   connmin = connSorted(end-300);
%   conn(conn<connmin) = -Inf;
%   connmax=max(conn(:));
%   
%   for k=1:66
%     
%     plot3(roi_center(k,1),roi_center(k,2),roi_center(k,3),'o','MarkerSize',13,'color',cm(:,k), 'MarkerFaceColor', cm(:,k),'LineWidth', 2)
%     for m=1:k-1
%       if 5*(conn(m,k)-connmin)/(connmax-connmin)>0
%         h=line([roi_center(k,1) roi_center(m,1)]',[roi_center(k,2) roi_center(m,2)]',[roi_center(k,3) roi_center(m,3)]','color','k','LineWidth',5*(conn(m,k)-connmin)/(connmax-connmin));
%       end
%     end
%   end
  hold off;
  
  
  clusterConnmat(logical(eye(size(clusterConnmat)))) = 0;
%   figure(4)
%   imagesc(clusterConnmat);
  

  
  
  for t=1:36
    az=az+10;
    view(az,el);
    drawnow;
    if doCreateVideo
      writeVideo(daObj,getframe(figh)); %use figure, since axis changes size based on view
    end
  end
end


if doCreateVideo
  close(daObj);
end