data = dataPaths();
path_results = fullfile(data.resultsdir, 'rgast/JansenRit/JR_Paper/DriverVarFreq1_inp_out');

results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

results.coh = squeeze(reshape(results.all_coh, paramSizes));

strengths = cell2mat(results.paramValues{3});
freqs = cell2mat(results.paramValues{4});


imagesc(squeeze(mean(results.coh,1)))
set(gca, 'XTick', 1:length(freqs), 'XTickLabel', freqs)
set(gca, 'YTick', 1:length(strengths), 'YTickLabel', strengths)