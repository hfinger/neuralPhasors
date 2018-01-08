
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

%% load freesurfer clustering:
tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/fsroibyvoxel.mat']);
clusterIdsFs = tmp.fsroi(useVoxelIdx);

clusterConnmatFs = load('/net/store/nbp/projects/phasesim/workdir/Holger/clusteredByFsROIs.mat');
clusterConnmatFs = clusterConnmatFs.clusterConnmat;
    
resortIds = load('/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIds.mat');
resortIds = resortIds.resortIds;
clusterIdsFsResorted = zeros(size(clusterIdsFs));
for k=1:66
  clusterIdsFsResorted(find(clusterIdsFs==resortIds(k))) = k;
end
clusterIdsFs = clusterIdsFsResorted;
clusterConnmatFs = clusterConnmatFs(resortIds,resortIds);

voxelIndexCoordsFs = zeros(66,3);
for k=1:66
  voxelIndexCoordsFs(k,:) = mean(usedVoxelCoord(find(clusterIdsFs==k),:),1);
end

%%
doReplaceByLowRes = false;
normByTargetROISize = true;
doCreateVideo = false;

%%
if doCreateVideo
  daObj=VideoWriter('connMat_basedOnFs_NormByProduct_climFix'); %for default video format. 
  daObj.FrameRate=6;
  open(daObj);
end

for i=66 %[66:200 202:2:400 404:4:600 608:8:800 816:16:1000 1000]
  
  %%
  if i==66
    voxelIndexCoords = voxelIndexCoordsFs;
    clusterIds = clusterIdsFs;
    clusterConnmat = clusterConnmatFs;
    
  else
    
    if i==67
      voxelIndexCoordsPrevious = voxelIndexCoordsFs;
    else
      tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fscos/clusterCenter' num2str(i-1) '.mat']);
      voxelIndexCoordsPrevious = tmp.clusterCenter;
    end
    
    tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fscos/clusterCenter' num2str(i) '.mat']);
    voxelIndexCoords = tmp.clusterCenter;
    tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fscos/fscosCluster' num2str(i) '.mat']);
    clusterIds = tmp.Cluster;
    tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/fscos/clusterConnmat' num2str(i) '.mat']);
    clusterConnmat = tmp.clusterConnmat;
    
  end
  
  
  %% resort to fs sorting:
  clusterIdsResorted = zeros(size(clusterIds));
  counter = 1;
  resortIdsCurrent = zeros(1,i);
  for k=1:66
    clusterIdsWithinFsRoiK = unique(clusterIds(clusterIdsFs==k));
    for j=1:length(clusterIdsWithinFsRoiK)
      resortIdsCurrent( counter) =  clusterIdsWithinFsRoiK(j);
      clusterIdsResorted(clusterIds==clusterIdsWithinFsRoiK(j)) = counter;
      counter = counter + 1;
    end
  end
  clusterIds = clusterIdsResorted;
  clusterConnmat = clusterConnmat(resortIdsCurrent,resortIdsCurrent);
  
  
  %%
  uniqueRoiIds = unique(clusterIds);
  numVoxInRoi = histc(clusterIds,uniqueRoiIds);
  clusterConnmat(logical(eye(length(uniqueRoiIds),length(uniqueRoiIds)))) = 0;
  oursProb = bsxfun(@rdivide, clusterConnmat, sum(clusterConnmat,2));
  if normByTargetROISize
    oursProb = bsxfun(@rdivide, oursProb, numVoxInRoi');
  end
  oursProb = oursProb + oursProb';
  
  %%
  
  if doReplaceByLowRes
    
    uniqueRoiIds = unique(clusterIdsFs);
    numVoxInRoi = histc(clusterIdsFs,uniqueRoiIds);
    clusterConnmatFs(logical(eye(length(uniqueRoiIds),length(uniqueRoiIds)))) = 0;
    oursProb = bsxfun(@rdivide, clusterConnmatFs, sum(clusterConnmatFs,2));
    if normByTargetROISize
      oursProb = bsxfun(@rdivide, oursProb, numVoxInRoi');
    end
    oursProb = oursProb + oursProb';
    
    fsId = zeros(size(clusterConnmat,1),1);
    for k=1:size(clusterConnmat,1)
      % search for the corresponding fs id
      fsId(k) = clusterIdsFs(find(clusterIds==k,1,'first'));
    end
    
    oursProbNew = zeros(size(clusterConnmat));
    for k=1:66
      ids = find(fsId==k);
      for j=1:66
        ids2 = find(fsId==j);
        oursProbNew(ids,ids2) = oursProb(k,j);
      end
    end
    
    oursProb = oursProbNew;
    
  end
  
  %%
  
  figh = figure(1);
  set(figh,'units','pixels','Position',[50 50 1000 1000]);
  imagesc(oursProb)
  cmap = colormap('jet');
%   set(gca,'clim',[0 1]);
  set(gca,'clim',[0 0.005]);
%   colormap(flipud(cmap))
  axis square;
  set(gca,'XTickLabel','')
  set(gca,'YTickLabel','')
  colorbar
  pause(0.01)
  title(['number ROI: ' num2str(i)]);
  drawnow;
  
  if doCreateVideo
    writeVideo(daObj,getframe(figh)); %use figure, since axis changes size based on view
  end
  
  
  %%
  if false%i>67
    
    
    indexOfSplitCluster = find(sum(abs(voxelIndexCoordsPrevious-voxelIndexCoords(1:end-1,:)),2));
    
    clusterOfVertex = clusterIds(voxIdOfVertex);
    clusterOfFace = clusterIds(voxIdOfFace);
    
    %   cdataPerVertex.cdata = cdataPerROI(clusterOfVertex,:);
    
    
    figure(3);
    clf;
    
    cdataPerVertex.cdata = repmat([0.5 0.5 0.5],size(g.vertices,1),1);
    %   hg=plot(g,cdataPerVertex);
    %   alpha(hg,0.1);
    %   shading flat
    %   view(-9.8705,72.0327);
    
    hold on;
    
    
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
    alphaVal = 1*ones(size(g.vertices,1),1);
    alphaVal([vertices1; vertices2]) = 1;
    shading flat
    set(hg,'FaceAlpha','interp')
    set(hg,'FaceVertexAlphaData',alphaVal)
    %   alpha(hg,0.5);
    
    
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
    
    plot3(clusterCenterCoord(:,1),clusterCenterCoord(:,2),clusterCenterCoord(:,3),'o','MarkerFaceColor',[0 0 0],'markersize',10)
    
    
    
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
    
    for t=1:30
      az=az+3;
      view(az,el);
      pause(0.05)
    end
    
  end
  
  
end


if doCreateVideo
  close(daObj);
end