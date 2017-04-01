
resultPath = '/net/store/nbp/projects/phasesim/results/rgast/JR_Driver_LTEffects';

files = dir(fullfile(resultPath,'JansenRitResults*'));
fnames = {files.name};
SID = 1;
data = load(strcat(resultPath, '/', fnames{SID}));

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

drivRange = 30; %data.simResult.sim.drivRange;
drivScale = 100; %data.simResult.sim.drivScale;
drivPos = data.simResult.sim.drivPos;
drivStrength = normpdf(ED(drivPos,:), 0, drivRange);
drivStrength = drivStrength - min(drivStrength);
drivConn = (drivStrength/max(drivStrength)) * drivScale;

resortRegions = 1;
plotColoredBrain(drivConn', resortRegions)