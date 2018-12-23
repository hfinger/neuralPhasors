data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['rgast/JansenRit/JR_Paper/' my_foldername]);

results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

results.coh = squeeze(reshape(results.all_coh, paramSizes));
results.FC = squeeze(reshape(results.all_FC, paramSizes));
results.corr_SimFC = squeeze(reshape(results.all_corr_SimFC, paramSizes));
results.corr_SimFC = cellfun(@(x) x{1}, results.corr_SimFC);

k = cell2mat(results.paramValues{1});
v = cell2mat(results.paramValues{2});

imagesc(results.corr_SimFC')
set(gca, 'XTick', 1:length(k), 'XTickLabel', k)
set(gca, 'YTick', 1:length(v), 'YTickLabel', v)
xlabel('k')
ylabel('v')