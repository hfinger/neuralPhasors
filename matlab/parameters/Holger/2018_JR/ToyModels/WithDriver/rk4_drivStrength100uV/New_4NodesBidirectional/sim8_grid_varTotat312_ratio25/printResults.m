clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ToyModels/WithDriver/rk4_drivStrength100uV/New_4NodesBidirectional/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

phase_offsets = [0:0.0625:1]*2*pi;
d12 = 0:10:200;

%% average over all sub sims:
all_FC = cell2mat(permute(cellfun(@(x) x{1},results.all_FC, 'UniformOutput', false),[1 3 2]));
all_FC = all_FC(33:end, 33:end, :);
num_nodes = 4;
num_repeat = 16;
FC_split = zeros(num_nodes, num_nodes, size(all_FC,3), num_repeat);
for k=1:num_repeat
  start_idx = (k-1)*num_nodes + 1;
  end_idx = k*num_nodes;
  FC_split(:,:,:,k) = all_FC(start_idx:end_idx, start_idx:end_idx, :);
end

FC_split = reshape(FC_split, [num_nodes, num_nodes, paramSizes, num_repeat]);

FC_split = squeeze(mean(FC_split, 5));
FC_split = permute(FC_split, [3, 4, 1, 2]);

%%
figure(1)
plotCohOfNodes(FC_split, 1, 2, d12, phase_offsets, path_results, 0, 0.5)
figure(2)
plotCohOfNodes(FC_split, 1, 3, d12, phase_offsets, path_results, 0, 0.5)
figure(3)
plotCohOfNodes(FC_split, 2, 4, d12, phase_offsets, path_results, 0, 0.5)
figure(4)
plotCohOfNodes(FC_split, 3, 4, d12, phase_offsets, path_results, 0, 0.5)
figure(5)
plotCohOfNodes(FC_split, 1, 4, d12, phase_offsets, path_results, 0, 0.5)

%%
PA1 = min( FC_split(:, :, 1,2), FC_split(:, :, 2,4));
PA2 = min( FC_split(:, :, 1,3), FC_split(:, :, 3,4));

%%
figure(6)
imagesc(d12,phase_offsets,PA1')
set(gca,'clim',[0; 0.3]);
set(gca,'TickLength',[0 0]);
set(gca,'xTick',[0, 100]);
set(gca,'yTick',[0]);
xlabel('distance [mm]')
ylabel('stimulation phase offset [rad]')
title('PA1')
set(gca,'YDir','normal')
colorbar;
saveas(gcf,fullfile(path_results,'PA1.png'))
saveas(gcf,fullfile(path_results,'PA1.pdf'))

%%
figure(7)
imagesc(d12,phase_offsets,PA2')
set(gca,'clim',[0; 0.3]);
set(gca,'TickLength',[0 0]);
set(gca,'xTick',[0, 100]);
set(gca,'yTick',[0]);
xlabel('distance [mm]')
ylabel('stimulation phase offset [rad]')
title('PA2')
set(gca,'YDir','normal')
colorbar;
saveas(gcf,fullfile(path_results,'PA2.png'))
saveas(gcf,fullfile(path_results,'PA2.pdf'))

%% calc IPSF:
[maxCoh, maxInd] = max(FC_split(:,1:16, 1, 4)');
[minCoh, minInd] = min(FC_split(:,1:16, 1, 4)');
stdCoh = std(FC_split(:,1:16, 1, 4)');
IPSF = maxCoh - minCoh;

figure(8);
set(gcf, 'Position',  [100, 100, 500, 250])
clf;
plot(d12,[IPSF]')
set(gca,'xTick',[0, 100]);
%set(gca,'yTick',[0.0, 0.1, 0.2]);
xlabel('distance [mm]')
ylabel('IPSF')
%set(gca,'ylim',[0.0; 0.2])
saveas(gcf,fullfile(path_results,'ipsf.png'))
saveas(gcf,fullfile(path_results,'ipsf.pdf'))