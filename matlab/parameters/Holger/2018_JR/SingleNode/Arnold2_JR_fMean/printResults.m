clear all;
data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/SingleNode/' my_foldername]);

results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

results.coh = squeeze(reshape(results.all_coh, paramSizes));

strengths = cell2mat(results.paramValues{2});
freqs = cell2mat(results.paramValues{3});

%%
figure(1);
imagesc(squeeze(results.coh))
set(gca, 'XTick', 1:2:length(freqs))
set(gca, 'XTickLabel', freqs(1:2:end))
set(gca, 'YTick', 1:length(strengths), 'YTickLabel', strengths)
colorbar;
xlabel('Freq [Hz]')
ylabel('Driver Strength [mV]')
title('coherence to driven node')
saveas(gcf,'coh_to_single_driven_node.png')