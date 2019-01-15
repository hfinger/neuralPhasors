clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ToyModels/WithDriver/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

%%
stim_pair_coh = zeros(length(results.all_FC),16);
for k=1:length(results.all_FC)
    FC = results.all_FC{k}{1}(33:end,33:end);
    FC_per_SPO = diag(FC,1);
    stim_pair_coh(k,:) = FC_per_SPO(1:2:31);
end

%%
imagesc(stim_pair_coh)