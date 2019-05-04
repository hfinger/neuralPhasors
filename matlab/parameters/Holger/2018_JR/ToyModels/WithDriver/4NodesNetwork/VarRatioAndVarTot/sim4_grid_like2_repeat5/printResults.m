clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ToyModels/WithDriver/4NodesNetwork/VarRatioAndVarTot/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

phase_offsets = [0:0.0625:1]*2*pi;
d213 = 0:10:200;
d24 = 0:10:100;

all_FC = cellfun(@(x) x{1}, results.all_FC, 'UniformOutput', false);
all_FC = permute(all_FC, [3, 4, 5, 1, 2]);
all_FC = cell2mat(all_FC);

all_FC = reshape(all_FC, [96, 96, length(d213), length(d24), 5]);
all_FC = mean(all_FC,5);
all_FC = reshape(all_FC, [96, 96, length(d213)*length(d24)]);

paramSizes = [length(d213), length(d24)];

%%
stim_pair_coh = zeros(size(all_FC,3),16);
PA1 = zeros(size(all_FC,3),16);
PA2 = zeros(size(all_FC,3),16);
for k=1:size(all_FC,3)
    FC = all_FC(33:end,33:end, k);
    diag1 = diag(FC,1);
    diag2 = diag(FC,2);
    diag3 = diag(FC,3);
    coh12 = diag1(1:4:end);
    coh23 = diag1(2:4:end);
    coh34 = diag1(3:4:end);
    coh13 = diag2(1:4:end);
    coh24 = diag2(2:4:end);
    coh14 = diag3(1:4:end);
    PA1(k,:) = min(coh12, coh13);
    PA2(k,:) = min(coh24, coh34);
    stim_pair_coh(k,:) = coh23;
end

[maxCoh, maxInd] = max(stim_pair_coh(:,1:16)');
[minCoh, minInd] = min(stim_pair_coh(:,1:16)');
stdCoh = std(stim_pair_coh(:,1:16)');
IPSF = maxCoh - minCoh;

PA1_at_max_SPO = diag(PA1(:,maxInd));
PA2_at_max_SPO = diag(PA2(:,maxInd));

stim_pair_coh_doubles = cat(2, stim_pair_coh, stim_pair_coh(:,1));
stim_pair_coh_doubles = reshape(stim_pair_coh_doubles, [paramSizes, 17]);
IPSF = reshape(IPSF, [paramSizes]);
PA1_at_max_SPO = reshape(PA1_at_max_SPO, [paramSizes]);
PA2_at_max_SPO = reshape(PA2_at_max_SPO, [paramSizes]);

PA1 = cat(2, PA1, PA1(:,1));
PA2 = cat(2, PA2, PA2(:,1));
PA1 = reshape(PA1, [paramSizes, 17]);
PA2 = reshape(PA2, [paramSizes, 17]);

diffPA = bsxfun(@rdivide, (PA1 - PA2), (PA1 + PA2));

fractionSPO_PA1_active = mean(diffPA(:, :, 1:16) > 0, 3);

%%
figure(1)

subplot(2,3,1)
imagesc(d24,d213,IPSF)
set(gca,'clim',[0; 0.6])
xlabel('node on P2 [mm]')
ylabel('total distance P1 [mm]')
title('IPSF')
set(gca,'YDir','normal')
colorbar;

%%
subplot(2,3,2)
imagesc(d24,d213,PA1_at_max_SPO)
set(gca,'clim',[0; 0.5])
xlabel('node on P2 [mm]')
ylabel('total distance P1 [mm]')
title('PA1 at max SPO')
set(gca,'YDir','normal')
colorbar;

%%
subplot(2,3,3)
imagesc(d24,d213,PA2_at_max_SPO)
set(gca,'clim',[0; 0.5])
xlabel('node on P2 [mm]')
ylabel('total distance P1 [mm]')
title('PA2 at max SPO')
set(gca,'YDir','normal')
colorbar;


%%
subplot(2,3,4)
imagesc(d24, phase_offsets, squeeze(diffPA(4,:,:))')
set(gca,'clim',[-0.5; 0.5])
xlabel('node on P2 [mm]')
ylabel('phase offset [rad]')
title('normalized difference between PA1 and PA2')
set(gca,'YDir','normal')
colorbar;

%%
subplot(2,3,5)
imagesc(d213, phase_offsets, squeeze(diffPA(:,9,:))')
set(gca,'clim',[-0.5; 0.5])
xlabel('total distance P1 [mm]')
ylabel('phase offset [rad]')
title('normalized difference between PA1 and PA2')
set(gca,'YDir','normal')
colorbar;

%%
subplot(2,3,6)
imagesc(d24,d213,fractionSPO_PA1_active)
set(gca,'clim',[0; 1])
xlabel('node on P2 [mm]')
ylabel('total distance P1 [mm]')
title('fraction of SPOs where PA1 > PA2')
set(gca,'YDir','normal')
colorbar;
