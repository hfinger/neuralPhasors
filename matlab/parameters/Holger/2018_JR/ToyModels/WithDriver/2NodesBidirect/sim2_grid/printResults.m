clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ToyModels/WithDriver/2NodesBidirect/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

phase_offsets = [0:0.0625:1]*2*pi;
d12 = 0:5:200; % distance between nodes at the ends of chain

%%
stim_pair_coh = zeros(length(results.all_FC),16);
for k=1:length(results.all_FC)
    FC = results.all_FC{k}{1}(33:end,33:end);
    FC_per_SPO = diag(FC,1);
    stim_pair_coh(k,:) = FC_per_SPO(1:2:31);
end
stim_pair_coh = cat(2, stim_pair_coh, stim_pair_coh(:,1));

%%
imagesc(phase_offsets,d12,stim_pair_coh)
xlabel('stimulation phase offset [rad]')
ylabel('distance [mm]')
set(gca,'YDir','normal')
colorbar;