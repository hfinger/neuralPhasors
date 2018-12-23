data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['rgast/JansenRit/JR_Paper/' my_foldername]);

results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

results.coh = squeeze(reshape(results.all_coh, paramSizes));

strengths = cell2mat(results.paramValues{2});
freqs = cell2mat(results.paramValues{3});


imagesc(squeeze(results.coh))
set(gca, 'XTick', 1:length(freqs), 'XTickLabel', freqs)
set(gca, 'YTick', 1:length(strengths), 'YTickLabel', strengths)
xlabel('Freq [Hz]')
ylabel('Driver Strength [mV]')