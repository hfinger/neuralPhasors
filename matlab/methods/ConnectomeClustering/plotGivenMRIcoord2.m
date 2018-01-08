
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
doCreateVideo = true;
condition = 'fullcon';
i=500;


%%
tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/' condition '/clusterCenter' num2str(i-1) '.mat']);
voxelIndexCoordsPrevious = tmp.clusterCenter;

tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/' condition '/clusterCenter' num2str(i) '.mat']);
voxelIndexCoords = tmp.clusterCenter;
if strcmp(condition,'fscon') || strcmp(condition,'fullcon')
  type = 'conn';
else
  type = 'cos';
end
if strcmp(condition(1:2),'fs')
  type = ['fs' type];
else
  type = ['full' type];
end
tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/' condition '/' type 'Cluster' num2str(i) '.mat']);
clusterIds = tmp.Cluster;
tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/Clusters/' condition '/clusterConnmat' num2str(i) '.mat']);
clusterConnmat = tmp.clusterConnmat;

indexOfSplitCluster = find(sum(abs(voxelIndexCoordsPrevious-voxelIndexCoords(1:end-1,:)),2));

clusterOfVertex = clusterIds(voxIdOfVertex);
clusterOfFace = clusterIds(voxIdOfFace);

%   cdataPerVertex.cdata = cdataPerROI(clusterOfVertex,:);





%%
if doCreateVideo
  daObj=VideoWriter(['topographyClusterVoxels1_iter' num2str(i) '_' condition]); %for default video format. 
  daObj.FrameRate=10;
  open(daObj);
end

colors = [0 0 1; 0 1 0];

%%

for k=1:20
  voxIds = find(clusterIds==k);
  voxCoordsCluster = usedVoxelCoord(voxIds,:);
  
  figh = figure(1);
  clf;
  hg = plot(g);
  alpha(hg,0.1);
  
  set(gca,'xtick',[],'ytick',[])
  
  hold on;
  hold on;
  plot3(voxCoordsCluster(:,1),voxCoordsCluster(:,2),voxCoordsCluster(:,3),'o','MarkerFaceColor',colors(mod(k,2)+1,:),'markersize',3)
  hold off;
  
  %set(gca,'xtick',[],'ytick',[])
  
  for t=1:36
    az=az+10;
    view(az,el);
%     pause(0.5)
    drawnow;
    if doCreateVideo
      writeVideo(daObj,getframe(figh)); %use figure, since axis changes size based on view
    end
  end
  
end

%%
if doCreateVideo
  close(daObj);
end