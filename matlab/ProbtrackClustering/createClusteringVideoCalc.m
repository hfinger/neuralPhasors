function [  ] = createClusteringVideoCalc( clustRange, CompSimPath, OutputPath, clusterIdPath )
%CREATECLUSTERINGVIDEOCALC Summary of this function goes here
%   Detailed explanation goes here

%% load tractography data:
path_prefix = '/net/store/nbp/projects/phasesim/databases/DTI_subject01/';
freesurfer_roi_ids = load([path_prefix 'tractographyData/freesurfer_roi_ids.mat']);
freesurfer_roi_ids = freesurfer_roi_ids.freesurfer_roi_ids;
header = cbiReadNiftiHeader([path_prefix 'surfaceFsData/compl_fs_mask.nii']);
g = gifti([path_prefix 'surfaceFsData/ca01_1_structcortex_8196.surf.gii']);

voxel_coords = load([CompSimPath 'FinalCoord/FinalCoord1.mat']);
voxel_coords = voxel_coords.FinalCoord;

%% create output folder
if ~exist(OutputPath, 'dir')
  mkdir(OutputPath);
end

ClusterIdsforVoxel = load(clusterIdPath);
largestClusterId = ClusterIdsforVoxel.largestClusterId;
ClusterIdsforVoxel = ClusterIdsforVoxel.stepWiseClusters;

figh = figure(1);
set(figh, 'Position', [100, 100, 330, 300]);
clf;

% load brain surface and plot it:
hg = plot(g);
alpha(hg,0); % is set to invisible
set(gca,'xtick',[],'ytick',[])
hold on;
% h1 = plotVoxelsInBrain(header, voxel_coords);

finalNumClust = max(clustRange);
finalClusterIdsForVoxel = ClusterIdsforVoxel(:,finalNumClust);
h = cell(finalNumClust,1);

%% create matrix that displays per iteration the cluster id per final clusters:
clusterSrcDst = zeros(finalNumClust,finalNumClust);
clusterSrcDst(end,:) = 1:finalNumClust;
for k=finalNumClust:-1:2
  lastRow = clusterSrcDst(k,:);
  lastRow(k) = largestClusterId(k);
  clusterSrcDst(k-1,:) = lastRow;
end
clusterSrcDst(1,:) = 1;
for k=2:finalNumClust
  row = clusterSrcDst(k,:);
  for i=k+1:finalNumClust
    row(i) = row(row(i));
  end
  clusterSrcDst(k,:) = row;
end

%% initialize once plot the patch of each final cluster:
clusterColors = distinguishable_colors(finalNumClust);
for k=1:finalNumClust
  h{k} = plotVoxelsInBrain(header, voxel_coords(finalClusterIdsForVoxel==k,:));
  set(h{k},'MarkerFaceColor',clusterColors(k,:),'markersize',4);
end

%% video
az=0;
maxElevation=40;
doCreateVideo = true;
if doCreateVideo
  daObj=VideoWriter(fullfile(OutputPath,'clustering_video')); %for default video format. 
  daObj.FrameRate=20;
  axis vis3d;
  camva('manual');
  camva(5.5);
  open(daObj);
end

initFramesPerCluster=13;

%% now only recolor the clusters accordingly at each iteration:
frameIter = 1;
for clusterCount = clustRange
  
  framesPerCluster = ceil(initFramesPerCluster*0.985.^clusterCount);
  
  for k=1:finalNumClust
    clustId = clusterSrcDst(clusterCount,k);
    set(h{k},'MarkerFaceColor',clusterColors(clustId,:))
  end
  
  for t=1:round(framesPerCluster)
    az = az + 2;
    el = maxElevation * cos( frameIter / 200 );
    view(az,el);
    drawnow;
    if doCreateVideo
      writeVideo(daObj,getframe(figh)); %use figure, since axis changes size based on view
    end
    
    if mod(frameIter, 100)==0
      disp(['frameIter =' num2str(frameIter)])
    end
    frameIter = frameIter + 1;
  end
end

if doCreateVideo
  close(daObj);
end


end

