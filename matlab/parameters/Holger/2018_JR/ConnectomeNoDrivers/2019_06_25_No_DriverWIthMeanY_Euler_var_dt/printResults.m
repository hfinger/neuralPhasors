clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
my_path = ['Holger/2018_JR/ConnectomeNoDrivers/' my_foldername];
path_results = fullfile(data.resultsdir, my_path);

if ~exist(path_results)
  mkdir(path_results)
end

jobDesc = load( fullfile(data.workdir, my_path ,'temp_Connectome','jobDesc.mat') );

fname = fullfile( data.workdir, my_path, 'Connectome1.mat');
tmp = load( fname );

%%
figure(1)
plot(tmp.simResult.spectrumFreq, log(tmp.simResult.spectrumPower{1}))
xlabel('Frequency [Hz]')
ylabel('log power')
saveas(gcf, fullfile( path_results, 'psd_of_driver.png'))

%%
spectrumPowerNetwork = cell2mat(tmp.simResult.spectrumPower(2:34)');
spectrumPowerNetworkMean = mean(spectrumPowerNetwork, 2);

%% calculate frequency with maximum power:
[maxVal,maxIdx] = max(spectrumPowerNetworkMean);
freqWithMaxPower = tmp.simResult.spectrumFreq(maxIdx);
disp(freqWithMaxPower)

%% plot spectrum of full network
figure(2)
plot(tmp.simResult.spectrumFreq, log(spectrumPowerNetworkMean))
xlabel('Frequency [Hz]')
ylabel('log power spectral density')
title(['freqWithMaxPower = ' num2str(freqWithMaxPower) ' Hz'])
saveas(gcf, fullfile( path_results, 'average_psd_of_all_nodes.png'))

%% plot spectrum of single nodes in network:
figure(3);
clf;
hold on;
for k=1:33
  smooth_filter_width = 50;
  smooth_filter = ones(1, smooth_filter_width)/smooth_filter_width;
  smoothed_psd = filter(smooth_filter, 1, log(spectrumPowerNetwork(:,k)));

  plot(tmp.simResult.spectrumFreq, (k-1)*2 + smoothed_psd)
end
xlabel('Frequency [Hz]')
ylabel('log power spectral density (with offsets for visibility)')
hold off;
saveas(gcf, fullfile( path_results, 'psd_of_each_node.png'))


%%

tmp = triu(ones(34,34), 1);
tmp(:,1) = 0;
tmp(1,:) = 0;
triuIdx = find(tmp(:));

results = cell(1,4);
for k=1:4
    results{k} = load(fullfile(data.workdir, my_path ,['/Connectome' num2str(k) '.mat']));
end


%%
coh_vals = cell(1,length(results));
for k=1:length(results)
  coh_vals{k} = results{k}.simResult.FC{1}(triuIdx);
  
  figure(k+10)
  imagesc(results{k}.simResult.FC{1})
end

coh_vals = cell2mat(coh_vals);
%%
all_corrs = corr(coh_vals);

figure(1);
imagesc(all_corrs)

%% calc mean abs differences:
meanAbsDiff = squeeze(mean(abs(bsxfun(@minus, coh_vals, permute(coh_vals,[1 3 2]))),1));

figure(1);
imagesc(meanAbsDiff)

%% compare timeseries:
node_idx = 5;
figure(2);
clf;
hold on;
for k=1:length(results)
  plot(k*1e-4 + results{k}.simResult.Y(node_idx,:))
  disp(std(results{k}.simResult.Y(node_idx,:)))
end
hold off;