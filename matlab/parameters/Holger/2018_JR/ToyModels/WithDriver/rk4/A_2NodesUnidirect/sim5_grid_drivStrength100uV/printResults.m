clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ToyModels/WithDriver/rk4/A_2NodesUnidirect/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

phase_offsets = [0:0.0625:1]*2*pi;
d12 = 0:5:200;

%%
stim_pair_coh = zeros(length(results.all_FC),16);
for k=1:length(results.all_FC)
    FC = results.all_FC{k}{1}(33:end,33:end);
    FC_per_SPO = diag(FC,1);
    stim_pair_coh(k,:) = FC_per_SPO(1:2:31);
end

stim_pair_coh = reshape(stim_pair_coh, [paramSizes, 16]);
stim_pair_coh = reshape(mean(stim_pair_coh,3), [paramSizes(2), 16]);
stim_pair_coh = cat(2, stim_pair_coh, stim_pair_coh(:,1));


%%
figure(1)
imagesc(d12,phase_offsets,stim_pair_coh')
%set(gca,'clim',[0.3; 0.7]);
set(gca,'TickLength',[0 0]);
set(gca,'xTick',[0, 200]);
set(gca,'yTick',[0]);
xlabel('distance [mm]')
ylabel('stimulation phase offset [rad]')
set(gca,'YDir','normal')
colorbar;
saveas(gcf,fullfile(path_results,'coh_vs_spo.png'))
saveas(gcf,fullfile(path_results,'coh_vs_spo.pdf'))

%% calc IPSF:
[maxCoh, maxInd] = max(stim_pair_coh(:,1:16)');
[minCoh, minInd] = min(stim_pair_coh(:,1:16)');
stdCoh = std(stim_pair_coh(:,1:16)');
IPSF = maxCoh - minCoh;

figure(2);
set(gcf, 'Position',  [100, 100, 500, 250])
clf;
plot(d12,[IPSF]')
set(gca,'xTick',[0, 200]);
set(gca,'yTick',[0, 0.2, 0.4, 0.6]);
xlabel('distance [mm]')
ylabel('IPSF')
set(gca,'ylim',[0; 0.6])
saveas(gcf,fullfile(path_results,'ipsf.png'))
saveas(gcf,fullfile(path_results,'ipsf.pdf'))
