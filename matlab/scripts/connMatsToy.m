
resultPath1 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_POvsScale';
resultPath2 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_ParamSearch';
resultPath3 = '/net/store/nbp/projects/phasesim/results/rgast';
resultPath4 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_PathTest';
resultPath6 = '/net/store/nbp/projects/phasesim/results/rgast/JR_Driver_LTEffects';

coherence_measure1 = 'drivPosCoh';
coherence_measure2 = 'Coherence';
coherence_measure3 = 'corr_SimFC';

clusters = 0;
dist_driver = 66;
useSP = 1;

%dataStruct = getArgsFromFiles(resultPath6,'JansenRitResults*',coherence_measure2,'drivPos','drivStart','drivDur','d','sampling','dt');
dataStruct = getArgsFromFiles(resultPath4,'JR_2Driver_PathTest*',{coherence_measure1, coherence_measure2},'p','drivPO','k','v','drivPos');

fnames = fieldnames(dataStruct);
dataTmp = dataStruct.(fnames{1});
%stimPos = dataTmp.drivPos;

%% plot k delay-weighted SWPs between driven regions
subject = 3;
data = dataStruct.(fnames{subject});

% set parameters
k = 3;
maxPathLength = 6;
p = data.p;
v = data.v;
targets = data.drivPos;
targets = targets - length(targets);
Coh = data.Coherence{1,1}(length(targets)+1:end,length(targets)+1:end);

% get delay-weighted SPWs
[paths, ~] = getDelayWeightedSWPs(targets, maxPathLength, p, v);

% compute pathlengths of SWPs based on Coherence matrix
pathlengths = zeros(1,k);
for i=1:k
    for j=1:length(paths{1,i})-1
        pathlengths(i) = pathlengths(i) + Coh(paths{1,i}(j),paths{1,i}(j+1));
    end
    pathlengths(i) = pathlengths(i) / (length(paths{1,i})-1);
end

% order the pathlengths and get indices
[pl_ordered,idx] = sort(pathlengths);
pl_ordered = fliplr(pl_ordered);
idx = fliplr(idx);

% order paths according to indices
for i=1:length(idx)
    paths_ordered{i} = paths{1,idx(i)};
end

% get node and edge values for path plotting
edgeVals = zeros(size(Coh));
for i=1:k
    for j=1:length(paths_ordered{1,i})-1
        edgeVals(paths_ordered{1,i}(j),paths_ordered{1,i}(j+1)) = Coh(paths_ordered{1,i}(j),paths_ordered{1,i}(j+1));
    end
end
nodeVals = zeros(1,size(Coh,1));
nodeVals(targets) = 1;

% plot nodes and edges in brain
edgeMin = 0.0001; % cut-off value for edgeVals (values below won't be plotted)
plotBrainConnectivity(nodeVals,edgeVals,edgeMin,1,1)

%% plot coherence between two driven regions over two dvs
logCoh = false;
clim = [0,1];
drivPosCohs = pltCohOver2DVs(dataStruct,{'drivPO','drivScale'},clim,logCoh);

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
