 function [  ] = createClusteringVideoFScalc( CompSimPath, OutputPath, FSClusterPath )
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

ClusterIdsforVoxel = load(fullfile(FSClusterPath, '01fsroiRemoveOverlap.mat'));
ClusterIdsforVoxel = ClusterIdsforVoxel.fsroi;

figh = figure(1);
set(figh, 'Position', [100, 100, 330, 300]);
clf;

% load brain surface and plot it:
hg = plot(g);
alpha(hg,0); % is set to invisible
set(gca,'xtick',[],'ytick',[])
hold on;

%% initialize once plot the patch of each final cluster:
clusterColors = distinguishable_colors(66);
for k=1:66
  h{k} = plotVoxelsInBrain(header, voxel_coords(ClusterIdsforVoxel==k,:));
  set(h{k},'MarkerFaceColor',clusterColors(k,:),'markersize',4);
end

%% video
az=0;
maxElevation=40;
doCreateVideo = true;
if doCreateVideo
  daObj=VideoWriter(fullfile(OutputPath,'fs_video')); %for default video format. 
  daObj.FrameRate=25;
  axis vis3d;
  camva('manual');
  camva(5.5);
  open(daObj);
end

%% now only recolor the clusters accordingly at each iteration:
frameIter = 1;
for t=1:720
  az = az + 2;
  el = maxElevation * cos( frameIter * pi / 360 );
  view(az,el);
  drawnow;
  if doCreateVideo
    writeVideo(daObj,getframe(figh)); %use figure, since axis changes size based on view
  end

  if mod(frameIter, 100)==0
    disp(['frameIter =' num2str(frameIter)])
  end
  frameIter = frameIter + 1;
  
  if frameIter == 360
    for k=1:66
      set(h{k},'markersize',1);
    end
  end
end

if doCreateVideo
  close(daObj);
end


end



