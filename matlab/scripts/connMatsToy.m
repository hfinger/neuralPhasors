
resultPath1 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/ToyModels/DelayEffects';
resultPath2 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_ParamSearch';
resultPath3 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_POvsScale/InpToExc'; %8_4, 11_8, 17_50, 21_32, 29_17, 30_7
resultPath4 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/2Driver_PathTest';
resultPath6 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/ToyModels/ArnoldTongue';
resultPath7 = '/net/store/nbp/projects/phasesim/results/rgast/JansenRit/NetworkParameters/InpToPyr';

measure1 = 'drivPosFC';
measure2 = 'FC';
measure3 = 'corr_SimFC';
measure4 = 'crossCorr';
measure5 = 'Y';
measure6 = 'freqs';

clusters = 0;
dist_driver = 66;
useSP = 1;

dataStruct = getArgsFromFiles(resultPath1,'JR_DelayEffects*',{measure2, measure6},'C','D','drivFreq','drivScale','drivPO','drivPos');

fnames = fieldnames(dataStruct);
dataTmp = dataStruct.(fnames{1});

%% build new data structures, dependent on entries in certain parameter

newStructs = newStruct(dataStruct,'C','D',[1,2]);

%% plot coherence between two driven regions over two dvs
dataStruct_tmp = newStructs.newStruct_3;
logCoh = false;
clim = [0,1];
pltType = 'im';
drivPos = dataTmp.drivPos;
targets = [drivPos];
cohIdx = 2;
baseline = [];
[Cohs, dvs] = pltCohOver2DVs(dataStruct_tmp,{'drivPO','D'},targets,1,clim,logCoh,pltType,cohIdx,baseline);

%% plot mutual coherence for driven regions over different target phase offsets
stimPos = dataTmp.drivPos;
POI = [2.35 5.]; % target phase offsets (does not need to match the exact simulation values)
mutualCoherence = cell(1,length(POI));
for p=1:length(POI)
    mutualCoherence{p} = plotMutualCoherence(dataStruct,POI(p),'drivScale',1);
end

%% plot difference in mutual coherence between max and min coherence PO in connectome

% calculate difference between mutual coherence at min and max PO
FCidx = 3;
[cohDiff, drivScales] = plotMutualCoherence(dataStruct,FCidx,'minmax','drivScale',0);

% load connectivity matrix
homotopeScaling = 0.1;
p = 1;
[C,~] = getConnectome(1,p,homotopeScaling,0);

% visualization parameters
m = max(max(cohDiff));
nodeRange = [0,m];
nodeScale = 20/m; % how strong should node size be scaled with nodeVal
nodeMinSize = 10; % minimum size of network nodes
CColor = 'k'; % edge color
CMin = 0.02; % cut-off value for edges
surfaceVisibility = 0.2; % transparency of brain surface (0 == no brain surface)

% for each target drivScale plot difference in mutual coherence over all nodes
targets = [0.6]; % target drivScales
for i=1:length(targets)
    idx = round(drivScales,3) == targets(i);
    figh = figure();
    nodeVals = [cohDiff(idx,:),zeros(1,33)];
    plotBrainConnectivity( nodeVals, C, nodeMinSize, nodeScale, nodeRange, CColor, CMin, [], surfaceVisibility, 1, 1 )
    set(figh, 'name',['Difference in mutual coherence for drivScale = ',num2str(targets(i))])
end

%% plot k delay-weighted SWPs between driven regions
data = dataStruct.(fnames{1});

% set parameters
k = 10;
kPlt = 1;
maxPathLength = 5;
targets = data.drivPos;
targets = targets - length(targets);

% load C
homotopeScaling = 0.1;
p = 1;
[C,~] = getConnectome(1,p,homotopeScaling,1);

threshold = 0.05;
invC = true;
[paths, ~] = getDelayWeightedSWPs(C, threshold, targets, maxPathLength, invC);

% target values for drivScale and phase offset
targetVals = [ 0.6 0.6 ; % target values for drivScale
                0. 1.7583]; % target values for drivPO

% visualization parameters
edgeMin = 0.001; % cut-off value for edges
nodeMin = 10; % minimum size of nodes
nodeScale = 20; % how strong node size should scale with nodeVal
invCoh = true; % invert coherence to calculate most active path or not
edgeScale = 10; % how strong edge size should scale with edge value
useNodeCoh = false; % whether to color nodes by their coherence with driven regions
pltEdgeVal = true; % if false, edge thickness will scale with activation of the whole path, instead of with edge activation
surfaceVisibility = 0.2; % determines transparency of the brain surface

% get node values
drivPOs = unique(targetVals(2,:));
drivScales = unique(targetVals(1,:));
nodeVals = zeros(length(drivScales),length(drivPOs),66);
FCidx = 3;
for i=1:length(drivPOs)
    [coherence, dv] = plotMutualCoherence(dataStruct,FCidx,drivPOs(i),'drivScale',0);
    for j=1:length(drivScales)
        nodeVals(j,i,1:33) = coherence(dv == round(drivScales(j),2),:);
    end
end

% plot paths in connectome
results = plotBrainInfoChannels(paths,dataStruct,FCidx,k,kPlt,{'drivScale','drivPO'},targetVals,nodeVals,nodeMin,nodeScale,edgeMin,edgeScale,pltEdgeVal,invCoh,useNodeCoh,surfaceVisibility);

% plot activation of most active path over phase offsets and driver scales
drivScales = unique(results.drivScale);
pathActivations = zeros(length(unique(results.drivPO)),length(drivScales));
yticklabels = cell(1,length(drivScales));
for d=1:length(drivScales)
    idx1 = results.drivScale == drivScales(d);
    [POs,idx2] = sort(results.drivPO(idx1));
    pl_tmp = results.pathActivations(idx1,1);
    pathActivations(:,d) = pl_tmp(idx2);
    yticklabels{d} = num2str(drivScales(d));
end
figure()
imagesc(pathActivations')
set(gca,'YTick',[1:length(drivScales)],'YTickLabel',yticklabels)
colorbar()
title('Activations of Shortest Paths')
xlabel('drivPO')
ylabel('drivScale')

%% Calculate mutual information for each pair of nodes
nWindows = 1;
nBins = 1000;
n = length(fnames);
MIs = cell(1,n);

for f=1:n
    data = dataStruct.(fnames{f});
    MIs{f} = mutualInformation(data,1,1,nWindows,nBins);
end

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
[bestMatch, bestMatchInd, corrs] = getBestTrial(dataStruct,'corr_SimFC');

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
