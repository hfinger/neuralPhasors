clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ConnectomeTwoDrivers/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

%%
stim_pair_coh = zeros(length(results.all_FC),1);
for k=1:length(results.all_FC)
    FC = results.all_FC{k}{1}(3:end, 3:end);
    stim_pairs = results.paramComb{2,k};
    stim_pair_coh(k) = FC(stim_pairs(1),stim_pairs(2));
end

%%
stim_pair_coh = squeeze(reshape(stim_pair_coh, paramSizes));

%%
[maxVals, maxInd] = max(stim_pair_coh, [], 2);
[minVals, minInd] = min(stim_pair_coh, [], 2);
hist(minInd)
IPSF = maxVals - minVals;
hist(IPSF);

%%
C = results.params.JansenRitConnectomePaper.C;
C = C' + C;
C = C > 0;

driver_node_idx = cell2mat(results.paramValues{2});
strengths = cell2mat(results.paramValues{3});
freqs = cell2mat(results.paramValues{4});

driver_coh = squeeze(mean(results.coh,1));
driver_coh_to_net = squeeze(mean(results.coh_driv_to_net,1));

%%
figure(1);
imagesc(driver_coh)
set(gca, 'XTick', 1:2:length(freqs))
set(gca, 'XTickLabel', freqs(1:2:end))
set(gca, 'YTick', 1:length(strengths), 'YTickLabel', strengths)
colorbar;
xlabel('Freq [Hz]')
ylabel('Driver Strength [mV]')
title('coherence to driven node')
saveas(gcf,fullfile( path_results, 'coh_to_driven_node.png'))

figure(6);
imagesc(driver_coh_to_net)
set(gca, 'XTick', 1:2:length(freqs))
set(gca, 'XTickLabel', freqs(1:2:end))
set(gca, 'YTick', 1:length(strengths), 'YTickLabel', strengths)
colorbar;
xlabel('Freq [Hz]')
ylabel('Driver Strength [mV]')
title('mean coherence to all network nodes')
saveas(gcf,fullfile( path_results, 'mean_coh_to_net.png'))

%%
target_freq = 11;
target_coh = 0.7;
driver_idx = 1;

target_freq_idx = find(freqs == target_freq);
target_driver_strength_idx = find(driver_coh(:,target_freq_idx) > target_coh, 1);

FC = results.FC{driver_idx, target_driver_strength_idx, target_freq_idx};
peakFreqsPerNode = results.freqs{driver_idx, target_driver_strength_idx, target_freq_idx};
cur_driver_node_idx = driver_node_idx(driver_idx);
coh_roi_with_driver = FC(2:end,1);
coh_roi_with_roi = FC(2:end,2:end);
coh_roi_with_net = (sum(coh_roi_with_roi,2)-1)/(size(coh_roi_with_roi,2)-1);

[B,I] = sort(coh_roi_with_driver);
coh_roi_with_driver = coh_roi_with_driver(I);
coh_roi_with_net = coh_roi_with_net(I);

figure(2)
plot([coh_roi_with_driver, coh_roi_with_net])

figure(3); 
imagesc(FC)

figure(4);
plot(peakFreqsPerNode)
title('peakFreqsPerNode');

%%
coh_driv_to_nn = zeros(size(results.FC));
for driver_idx=1:size(results.FC,1)
    cur_driver_node_idx = driver_node_idx(driver_idx);
    all_nn_idxs = find(C(cur_driver_node_idx,:));
    tmp = cellfun(@(x) x(1,2:end), results.FC(driver_idx, :, :), 'UniformOutput', false);
    coh_driv_to_nn(driver_idx, :, :) = cell2mat(cellfun(@(x) mean(x(all_nn_idxs)), tmp, 'UniformOutput', false));
end

figure(9);
imagesc(squeeze(mean(coh_driv_to_nn,1)))
set(gca, 'XTick', 1:2:length(freqs))
set(gca, 'XTickLabel', freqs(1:2:end))
set(gca, 'YTick', 1:length(strengths), 'YTickLabel', strengths)
colorbar;
xlabel('Freq [Hz]')
ylabel('Driver Strength [mV]')
title('mean coherence to nearest neigbors in network')
saveas(gcf,fullfile( path_results, 'mean_coh_to_nn.png'))
