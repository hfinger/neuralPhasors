clear all;
triuIds = find(triu(ones(66,66),1));

%% load our roiSize:
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
tmp = load(['/net/store/nbp/projects/phasesim/workdir/Arushi/20150615RecursiveNcut/fsroibyvoxel.mat']);
clusterIdsFs = tmp.fsroi(useVoxelIdx);
uniqueRoiIds = unique(clusterIdsFs);
numVoxInRoi = histc(clusterIdsFs,uniqueRoiIds);

%% load their tracking results:
theirsLoaded = load('/net/store/nbp/projects/phasesim/databases/SC_Bastian/dti_20141209_preprocessed.mat');
theirsRaw = theirsLoaded.ci{1};
roisize = theirsLoaded.roisize{1};

%remove NaN
theirsRaw(isnan(theirsRaw)) = 0;

%make symmetric
theirsSym = triu(theirsRaw,1)+triu(theirsRaw,1)';

%only the triangular part
theirsVec = theirsSym(triuIds);

%normalize for SAR
theirsNorm = bsxfun(@rdivide, theirsSym, sum(theirsSym,2));


%% load our tracking results:


oursRaw = load('/net/store/nbp/projects/phasesim/workdir/Holger/clusteredByFsROIs.mat');
oursRaw = oursRaw.clusterConnmat;

%remove diagonal entries:
oursRaw(logical(eye(66,66))) = 0;

%normalize 
% oursProb = bsxfun(@rdivide, oursRaw, sum(oursRaw,2));
oursProb = oursRaw;

% oursCorrected = bsxfun(@rdivide, oursProb, roisize);
oursCorrected = bsxfun(@rdivide, oursProb, numVoxInRoi');
% oursCorrected = bsxfun(@rdivide, oursProb, roisize * roisize');
% oursCorrected = oursProb;

%make symmetric
oursSym = oursCorrected + oursCorrected';

%correct for region sizes:
% roisize = theirsLoaded.roisize{1};
% oursSym = bsxfun(@rdivide, oursSym, (roisize*roisize'));

%only the triangular part
oursVec = oursSym(triuIds);

%normalize for SAR
oursNorm = bsxfun(@rdivide, oursSym, sum(oursSym,2));


%% correlate our and their SC:
c=corr(oursVec,theirsVec);
disp(['corr(theirs,ours)=' num2str(c)])

%% load FC:
load('/net/store/nbp/projects/phasesim/databases/SC_Bastian/eeg_20150114_controls_fs_funconn_lcmv_bponetrial_hilbert_3_30.mat')
FC=squeeze(mean(mean(coh_all(1,1:2,5:6,:,:,7),2),3));

%% calculate SAR model:
k=0.4:0.01:0.95; 
for i=1:length(k); 
  [oursCov, oursFC] = sar(oursNorm, k(i)); 
  oursPerf(i)=corr( FC(triuIds), oursFC(triuIds)); 
  
  [theirsCov, theirsFC] = sar(theirsNorm, k(i)); 
  theirsPerf(i)=corr( FC(triuIds), theirsFC(triuIds)); 
end;

disp(['our   performance SAR: ' num2str(max(oursPerf))])
disp(['their performance SAR: ' num2str(max(theirsPerf))])