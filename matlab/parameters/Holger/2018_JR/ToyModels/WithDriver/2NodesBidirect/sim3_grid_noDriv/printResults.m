clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ToyModels/WithDriver/2NodesBidirect/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

nbins = 32;
phase_offsets_hist_edges = linspace(0, 2*pi, nbins+1);
phase_offsets_hist_centers = phase_offsets_hist(1:end-1) + (phase_offsets_hist(2) - phase_offsets_hist(1)) / 2;
d12 = 0:5:200; % distance between nodes at the ends of chain

%%
stim_pair_coh = zeros(length(results.all_FC),nbins);
for k=1:length(results.all_FC)
    PhaseHist = results.all_FC{k}{2}(33:end,33:end,:);
    tmp = num2cell(PhaseHist,[1, 2]);
    tmp2 = cellfun(@(x) diag(x,1), tmp, 'UniformOutput', false);
    tmp3 = squeeze(cell2mat(tmp2));
    tmp4 = tmp3(1:2:31, :);
    tmp5 = sum(tmp4,1);
    tmp6 = tmp5 / sum(tmp5);
    stim_pair_coh(k,:) = tmp6;
end
stim_pair_coh = cat(2, stim_pair_coh, stim_pair_coh(:,1));

%%
imagesc(phase_offsets_hist_centers,d12,stim_pair_coh)
xlabel('phase offset [rad]')
ylabel('distance [mm]')
set(gca,'YDir','normal')
colorbar;