clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ConnectomeTwoDrivers/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

drivStrength = cell2mat(results.paramValues{3});
phaseOffset = cell2mat(results.paramValues{4});

%%
stim_pair_coh = zeros(lengthdifferences in Figure 6A ar(results.all_FC),1);
for k=1:length(results.all_FC)
    FC = results.all_FC{k}{1}(3:end, 3:end);
    stim_pairs = results.paramComb{2,k};
    stim_pair_coh(k) = FC(stim_pairs(1),stim_pairs(2));
end

%%
stim_pair_coh = squeeze(reshape(stim_pair_coh, paramSizes));
stim_pair_coh = mean(stim_pair_coh, 4);

%%
for k=1:6
    figure(k)
    imagesc(squeeze(stim_pair_coh(k,:,:)))
    set(gca,'YDir','normal')
    set(gca, 'XTick', 1:2:length(phaseOffset))
    set(gca, 'XTickLabel', phaseOffset(1:2:end))
    set(gca, 'YTick', 1:length(drivStrength), 'YTickLabel', drivStrength)
    set(gca, 'clim', [0, 1])
    colorbar;
    xlabel('Phase Offset')
    ylabel('Driver Strength [mV]')
    title('coherence between the two driven nodes')
    saveas(gcf,fullfile( path_results, ['coh_between_pair_' num2str(k) '.png']))
end

%%
[maxVals, maxInd] = max(stim_pair_coh, [], 3);
[minVals, minInd] = min(stim_pair_coh, [], 3);
hist(minInd)
IPSF = maxVals - minVals;
hist(IPSF);



