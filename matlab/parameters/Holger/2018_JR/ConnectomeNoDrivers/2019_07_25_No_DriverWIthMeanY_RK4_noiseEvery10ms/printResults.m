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
% 
% %%
% figure(1)
% plot(tmp.simResult.spectrumFreq, log(tmp.simResult.spectrumPower{1}))
% xlabel('Frequency [Hz]')
% ylabel('log power')
% saveas(gcf, fullfile( path_results, 'psd_of_driver.png'))
% 
% %%
% spectrumPowerNetwork = cell2mat(tmp.simResult.spectrumPower(2:34)');
% spectrumPowerNetworkMean = mean(spectrumPowerNetwork, 2);
% 
% %% calculate frequency with maximum power:
% [maxVal,maxIdx] = max(spectrumPowerNetworkMean);
% freqWithMaxPower = tmp.simResult.spectrumFreq(maxIdx);
% disp(freqWithMaxPower)
% 
% %% plot spectrum of full network
% figure(2)
% plot(tmp.simResult.spectrumFreq, log(spectrumPowerNetworkMean))
% xlabel('Frequency [Hz]')
% ylabel('log power spectral density')
% title(['freqWithMaxPower = ' num2str(freqWithMaxPower) ' Hz'])
% saveas(gcf, fullfile( path_results, 'average_psd_of_all_nodes.png'))
% 
% %% plot spectrum of single nodes in network:
% figure(3);
% clf;
% hold on;
% for k=1:33
%   smooth_filter_width = 50;
%   smooth_filter = ones(1, smooth_filter_width)/smooth_filter_width;
%   smoothed_psd = filter(smooth_filter, 1, log(spectrumPowerNetwork(:,k)));
% 
%   plot(tmp.simResult.spectrumFreq, (k-1)*2 + smoothed_psd)
% end
% xlabel('Frequency [Hz]')
% ylabel('log power spectral density (with offsets for visibility)')
% hold off;
% saveas(gcf, fullfile( path_results, 'psd_of_each_node.png'))


%%

tmp = triu(ones(34,34), 1);
tmp(:,1) = 0;
tmp(1,:) = 0;
triuIdx = find(tmp(:));

results = cell(1,7);
all_dt = zeros(1,length(results));
for k=1:length(results)
    results{k} = load(fullfile(data.workdir, my_path ,['/Connectome' num2str(k) '.mat']));
    all_dt(k) = results{k}.simResult.sim.dt;
end


%%
coh_vals = cell(1,length(results));
for k=1:length(results)
  coh_vals{k} = results{k}.simResult.FC{1}(triuIdx);
  
%   figure(k+10)
%   imagesc(results{k}.simResult.FC{1})
end

coh_vals = cell2mat(coh_vals);

%%
node_idx = 5;
clear timeseries;
for k=1:7
  timeseries(k,:) = results{k}.simResult.Y(node_idx,1:end);
end
all_corrs_raw = corr(timeseries');

figure(13);
imagesc(all_corrs_raw)
set(gca, 'XTick', 1:2:length(results))
set(gca, 'XTickLabel', all_dt(1:2:length(results)))
set(gca, 'YTick', 1:length(results))
set(gca, 'YTickLabel', all_dt)
set(gca,'YDir','normal')
xlabel('dt [sec]')
ylabel('dt [sec]')
colorbar;
saveas(gcf, fullfile( path_results, 'corr_between_timeseries.png'))

%%
all_corrs = corr(coh_vals);

figure(1);
imagesc(all_corrs)
set(gca, 'XTick', 1:2:length(results))
set(gca, 'XTickLabel', all_dt(1:2:length(results)))
set(gca, 'YTick', 1:length(results))
set(gca, 'YTickLabel', all_dt)
set(gca,'YDir','normal')
xlabel('dt [sec]')
ylabel('dt [sec]')
colorbar;
saveas(gcf, fullfile( path_results, 'corr_between_FCs.png'))

%% calc mean abs differences:
meanAbsDiff = squeeze(mean(abs(bsxfun(@minus, coh_vals, permute(coh_vals,[1 3 2]))),1));

figure(2);
imagesc(meanAbsDiff)
set(gca, 'XTick', 1:2:length(results))
set(gca, 'XTickLabel', all_dt(1:2:length(results)))
set(gca, 'YTick', 1:length(results))
set(gca, 'YTickLabel', all_dt)
set(gca,'YDir','normal')
xlabel('dt [sec]')
ylabel('dt [sec]')
colorbar;
saveas(gcf, fullfile( path_results, 'distance_between_FCs.png'))

%% suplementary figure 1:
node_idx = 5;
fig = figure(25);
clf;
set(gcf, 'Position',  [100, 100, 1500, 350])
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 0.2]);
%fig.PaperUnits = 'centimeters';
%fig.PaperPosition = [0 0 15.0 3.5];

subplot(1,3,1)
hold on;
times = (0:(length(results{1}.simResult.Y(node_idx,:))-1))/100;
for k=1:length(results)
  plot(times, k*1e-4 + results{k}.simResult.Y(node_idx,:))
end
set(gca,'Ytick',(1:length(results))*1e-4+1e-4)
set(gca,'Yticklabels',all_dt)
xlim([0, 5])
ylabel('dt [sec]')
xlabel('time [sec]')
hold off;

subplot(1,3,2)
imagesc(all_corrs_raw)
axis square;
set(gca, 'XTick', 1:2:length(results))
set(gca, 'XTickLabel', all_dt(1:2:length(results)))
set(gca, 'YTick', 1:length(results))
set(gca, 'YTickLabel', all_dt)
set(gca,'YDir','normal')
xlabel('dt [sec]')
ylabel('dt [sec]')
ax = gca;
ax.CLim = [0.85, 1];
colorbar;

subplot(1,3,3)
imagesc(all_corrs)
axis square;
set(gca, 'XTick', 1:2:length(results))
set(gca, 'XTickLabel', all_dt(1:2:length(results)))
set(gca, 'YTick', 1:length(results))
set(gca, 'YTickLabel', all_dt)
set(gca,'YDir','normal')
xlabel('dt [sec]')
ylabel('dt [sec]')
ax = gca;
ax.CLim = [0.55, 1];
colorbar;

saveas(gcf, fullfile( path_results, 'sup1.png'))
saveas(gcf, fullfile( path_results, 'sup1.pdf'))
