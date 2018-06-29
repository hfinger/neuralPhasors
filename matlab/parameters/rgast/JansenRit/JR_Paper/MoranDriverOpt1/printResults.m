data = dataPaths();
path_results = fullfile(data.resultsdir, 'rgast/JansenRit/JR_Paper/MoranDriverOpt1');

results = load(fullfile( path_results, 'all_coh.mat'), 'all_coh', 'paramComb', 'variableParams');

results.coh = (reshape(results.all_coh, [5, 8, 13]));

imagesc(squeeze(mean(results.coh,1)))