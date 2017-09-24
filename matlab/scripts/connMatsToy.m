
resultPath1 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_POvsScale';
resultPath2 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_ParamSearch';
resultPath3 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_POvsScale/nodes17_50';
resultPath4 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_PathTest';
resultPath6 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_InfoChannels';

coherence_measure1 = 'drivPosCoh';
coherence_measure2 = 'Coherence';
coherence_measure3 = 'corr_SimFC';
coherence_measure4 = 'crossCorr';

clusters = 0;
dist_driver = 66;
useSP = 1;

%dataStruct = getArgsFromFiles(resultPath6,'JansenRitResults*',coherence_measure2,'drivPos','drivStart','drivDur','d','sampling','dt');
dataStruct = getArgsFromFiles(resultPath3,'JR_2Driver_POvsScale_1750*',{coherence_measure2, coherence_measure1},'p','drivPO','k','v','drivPos','drivScale');

fnames = fieldnames(dataStruct);
dataTmp = dataStruct.(fnames{1});

%% plot mutual coherence for driven regions over different target phase offsets
stimPos = dataTmp.drivPos;
POI = [0,2];
mutualCoherence = cell(1,length(POI));
for p=1:length(POI)
    mutualCoherence{p} = plotMutualCoherence(dataStruct,POI(p),'drivScale',1);
end

%% plot difference in mutual coherence between max and min coherence PO in connectome

addpath('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e')

% get coherence between driven regions over POs and drivScales
drivPosCohs = pltCohOver2DVs(dataStruct,{'drivPO','drivScale'},0);

% calculate difference between mutual coherence at min and max PO
cohDiff = plotMutualCoherence(dataStruct,'minmax','drivScale',0,drivPosCohs');

% load connectivity matrix
path_SCmat = '/net/store/nbp/projects/phasesim/databases/avg_SC.mat';
path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
load(path_SCmat);
load(path_ResortIDs);
resortIds = [resortIdsMarlene, resortIdsMarlene + 33];
clear D

% reorder and rescale connectivity matrix
C = C(resortIds,:);
C = C(:,resortIds);
C = C + 0.1 * diag(ones(size(C,1)/2,1),size(C,1)/2) + 0.1 * diag(ones(size(C,1)/2,1),-size(C,1)/2);
C = bsxfun(@rdivide,C,sum(C.^2,2).^(1/2));

% for each drivScale plot difference in mutual coherence over all nodes
nodeScale = 30/max(max(cohDiff));
nodeMinSize = 10;
CColor = 'k';
CMin = 0.1;
surfaceVisibility = 0.1;

drivScales = [0:0.001:0.025];
for i=2:3:20
    m = max(cohDiff(i,:));
    nodeRange = [-m,m];
    figh = figure();
    plotBrainConnectivity( cohDiff(i,:), C, nodeMinSize, nodeScale, nodeRange, CColor, CMin, [], surfaceVisibility, 1, 1 )
    set(figh, 'name',['Difference in mutual coherence for drivScale = ',num2str(drivScales(i))])
end

%% plot k delay-weighted SWPs between driven regions
data = dataStruct.(fnames{1});

% set parameters
k = 10;
kPlt = 3;
maxPathLength = 4;
p = data.p;
v = data.v;
targets = data.drivPos;
targets = targets - length(targets);

% get delay-weighted SPWs
[paths, ~] = getDelayWeightedSWPs(targets, maxPathLength, p, v);

% get node values
targetVals = [ 0.001 0.007 0.013;0. 0. 0.2];
drivPosCohs = pltCohOver2DVs(dataStruct,{'drivScale','drivPO'},0);
mutualCoh1 = plotMutualCoherence(dataStruct,targetVals(2,1),'drivScale',0,drivPosCohs');
mutualCoh2 = plotMutualCoherence(dataStruct,targetVals(2,2),'drivScale',0,drivPosCohs');
nodeVals = zeros(size(mutualCoh1,1),2,size(mutualCoh1,2));
nodeVals(:,1,:) = mutualCoh1;
nodeVals(:,2,:) = mutualCoh2;
nodeVals = nodeVals(2:end,:,:);

edgeMin = 0.0005;
nodeMin = 10;
nodeScale = 20;
invCoh = 1;
edgeScale = 4;
results = plotBrainInfoChannels(paths,dataStruct,k,kPlt,{'drivScale','drivPO'},targetVals,nodeVals,nodeMin,nodeScale,edgeMin,edgeScale,invCoh);

drivScales = unique(results.drivScale);
pathActivations = zeros(length(unique(results.drivPO)),length(drivScales));
for d=1:length(drivScales)
    idx1 = results.drivScale == drivScales(d);
    [POs,idx2] = sort(results.drivPO(idx1));
    pl_tmp = results.pathActivations(idx1,1);
    pathActivations(:,d) = pl_tmp(idx2);
end
figure()
imagesc(pathActivations')
colorbar()
title('Activations of Shortest Paths')
xlabel('drivPO')
ylabel('drivScale')

%% plot coherence between two driven regions over two dvs
logCoh = false;
clim = [0,1];
meshplt = true;
drivPosCohs = pltCohOver2DVs(dataStruct,{'drivPO','drivScale'},clim,logCoh,meshplt);

%% plot coherence for different driver phase offsets and extract parameter set with highest variance in coherence over POs
logCoh = false;
clim = [0,1];
drivPosCohs = pltCohOverPhaseOffset(dataStruct,'drivScale',clim,logCoh);

drivPosCohSD = std(drivPosCohs,[],2);
[maxSD, idx] = max(drivPosCohSD);
highestSDFname = fnames{idx};
highestSDData = dataStruct.(fnames{idx});


%% plot phase difference between 2 target nodes
env.t_rm = 1; % time to start bandpassing the signal at [s]
env.sigBandpass(1).Fst1 = 25.5; % end stop band [Hz] 
env.sigBandpass(1).Fp1 = 28; % start pass band [Hz]
env.sigBandpass(1).Fp2 = 32; % [Hz] end pass band
env.sigBandpass(1).Fst2 = 34.5; % [Hz] start stop band
env.sigBandpass(1).Ast1 = 30; % frequency attenuation in first stopband
env.sigBandpass(1).Ap = 1; % passband ripples
env.sigBandpass(1).Ast2 = 30; % frequency attenuation in second stopband

drivEnd = dataTmp.drivStart + dataTmp.drivDur - dataTmp.d;
targets = [1, stimPos(1)];
pltWindow = [drivEnd - 50.5, drivEnd + 50.5];
cutOff = 1;

phaseDiff = plotPhaseDiffOverTime(dataStruct, targets, pltWindow, drivEnd, env);


%% look for best match between simulated and empirical FC
[bestMatch, bestMatchInd, corrs] = getBestTrial(dataStruct,'Coherence');

%% plot all connectivity matrices
results1 = plotConnMats(dataStruct, ['stimRange', 'stimScale'], coherence_measure2);

%% plot coherence over distance for each driver
for i=1:length(stimPos)
    [Coh_overDist] = plotCohOverDist(dataStruct, i, dist_driver, coherence_measure2, ['stimScale', 'stimRange'], useSP);
end

%% plot arnold tongue for each driver
for i=1:length(stimPos)
    [S] = arnTongue(dataStruct, [i,stimPos(i)], coherence_measure2, 1, {'drivPos', 'drivScale'}, clusters);
end

%% Average coherence over all simulations and plot mean coherence over chosen time windows
driverID = 2;
nDrivers = size(dataStruct.(fnames{1}).Coherence, 1);
Coherence_tmp = dataStruct.(fnames{1}).Coherence(driverID,:,:);
Coherence_tmp = squeeze(Coherence_tmp)';
CoherenceComplete = zeros(length(fnames),size(Coherence_tmp, 1), size(Coherence_tmp, 2));
for i = 1:length(fnames)
    Coherence_tmp = dataStruct.(fnames{1}).Coherence(driverID,:,:);
    Coherence_tmp = squeeze(Coherence_tmp)';
    CoherenceComplete(i,:,:) = Coherence_tmp;
end
CoherenceFinal = squeeze(mean(CoherenceComplete,1));
imagesc(CoherenceFinal)
colormap('jet')
colorbar()

TOI_Start = [1,21,101];
TOI_End = [20,100,length(CoherenceFinal)];
for j=1:length(TOI_Start)
    Coherence_Plt = squeeze(mean(CoherenceFinal(TOI_Start(j):TOI_End(j),:), 1));
    plotColoredBrain(Coherence_Plt(nDrivers+1:end)', 1)
end

%% show coherence and driven regions of single subject
SID = 1;
driverID = 1;
TOI_Start = 20;
TOI_End = 100;

path_ED = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/Stuct_connectivity_for_Holger.mat';
load(path_ED, 'roi_positions','struct_labels_correspnding_pos')
roi_coordinates = roi_positions(cell2mat(struct_labels_correspnding_pos),:);
path_ResortIDs = '/net/store/nbp/projects/phasesim/databases/SC_Bastian/resortIdsMarlene.mat';
load(path_ResortIDs);
resortIds = [resortIdsMarlene, resortIdsMarlene + 33];
roi_coordinates = roi_coordinates(resortIds,:);
electrodeCoordinatesTable = readtable('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/ca_electrodeLocations/ca03_EEGLocations.txt','Delimiter',' ','ReadVariableNames',false);
electrode_coordinates = table2array(electrodeCoordinatesTable(:,2:end));
ED = pdist2(electrode_coordinates, roi_coordinates);

SubjData = dataStruct.(fnames{SID});
Coherence_tmp = SubjData.Coherence;
Coherence_tmp = squeeze(Coherence_tmp)';
if size(Coherence_tmp, 1) > 1
    Coherence_tmp = mean(Coherence_tmp(TOI_Start:TOI_End,:), 1);
end
addpath('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e')
plotColoredBrain(Coherence_tmp(2:end)',1)

drivPos = SubjData.stimPos(driverID);
[~, electrodeIdx] = min(ED(:,drivPos));
drivStrength = normpdf(ED(electrodeIdx,:), 0, SubjData.stimRange);
drivStrength = drivStrength - min(drivStrength);
drivConn = (drivStrength/max(drivStrength)) * SubjData.stimScale;
plotColoredBrain(drivConn', 1)
