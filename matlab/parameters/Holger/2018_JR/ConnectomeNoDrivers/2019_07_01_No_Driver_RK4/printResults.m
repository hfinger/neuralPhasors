clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
my_path = ['Holger/2018_JR/ConnectomeNoDrivers/' my_foldername];
path_results = fullfile(data.resultsdir, my_path);

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
plot(tmp.simResult.spectrumFreq, 10*log10(spectrumPowerNetworkMean))
xlabel('Frequency [Hz]')
ylabel('log power spectral density')
title(['freqWithMaxPower = ' num2str(freqWithMaxPower) ' Hz'])
saveas(gcf, fullfile( path_results, 'average_psd_of_all_nodes.png'))

%% plot spectrum of single nodes in network:
fig3 = figure(3);
fig3.Renderer='Painters';

clf;
hold on;
for k=1:33
  smooth_filter_width = 50;
  smooth_filter = ones(1, smooth_filter_width)/smooth_filter_width;
  smoothed_psd = filter(smooth_filter, 1, 10*log10(spectrumPowerNetwork(:,k)));

  plot(tmp.simResult.spectrumFreq, (k-1)*2 + smoothed_psd)
  
  
  [~,maxIdx] = max(spectrumPowerNetwork(:,k));
  freqWithMaxPowerPerNode(k) = tmp.simResult.spectrumFreq(maxIdx);
  
end
xlabel('Frequency [Hz]')
ylabel('log power spectral density (with offsets for visibility)')
hold off;
saveas(gcf, fullfile( path_results, 'psd_of_each_node.png'))
saveas(gcf, fullfile( path_results, 'psd_of_each_node.pdf'))


%% plot spectrum of single nodes in network:
fig3 = figure(3);
fig3.Renderer='Painters';
cmap = lines;
clf;
hold on;
plot_nodes = [5, 10, 15];
for j=1:length(plot_nodes)
  k = plot_nodes(j);
  
  smooth_filter_width = 10;
  smooth_filter = ones(1, smooth_filter_width)/smooth_filter_width;
  smoothed_psd = filter(smooth_filter, 1, 10*log10(spectrumPowerNetwork(:,k)));

  plot_resolution = 10;
  subplot(length(plot_nodes),1,j)
  plot(tmp.simResult.spectrumFreq(1:plot_resolution:end), smoothed_psd((1:plot_resolution:end)), 'Color', cmap(j,:));
  if j==2
    ylabel('Power/frequency [db/Hz]')
  end
end
xlabel('Frequency [Hz]')
hold off;
saveas(gcf, fullfile( path_results, 'psd_of_3_noded.png'))
saveas(gcf, fullfile( path_results, 'psd_of_3_noded.pdf'))
