
resultPath = '/net/store/nbp/projects/phasesim/results/rgast/JR_Driver_LTEffects';

% load exemplary file
files = dir(fullfile(resultPath,'JansenRitResults*'));
fnames = {files.name};
SID = 1;
data = load(strcat(resultPath, '/', fnames{SID}));

% load roi positions and electrode positions and calculate eucl. distances
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

% create gaussian driver around 
drivRange = 30; %data.simResult.sim.drivRange;
drivScale = 100; %data.simResult.sim.drivScale;
drivPos = data.simResult.sim.drivPos;
[~, electrodeIdx] = min(ED(:,drivPos));
drivStrength = normpdf(ED(electrodeIdx,:), 0, drivRange);
drivStrength = drivStrength - min(drivStrength);
drivConn = (drivStrength/max(drivStrength)) * drivScale;

% plot driver around 
resortRegions = 1;
transparency = 1;
plotColoredBrain(drivConn', resortRegions, transparency)