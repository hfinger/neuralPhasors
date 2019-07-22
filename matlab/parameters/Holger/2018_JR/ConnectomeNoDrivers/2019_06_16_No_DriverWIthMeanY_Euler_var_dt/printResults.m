
tmp = triu(ones(34,34), 1);
tmp(:,1) = 0;
tmp(1,:) = 0;
triuIdx = find(tmp(:));

result_euler = cell(1,3);
for k=1:3
    result_euler{k} = load(['/home/hofinger/dev/matlab/workdir/Holger/2018_JR/ConnectomeNoDrivers/2019_06_16_No_DriverWIthMeanY_Euler_var_dt/Connectome' num2str(k) '.mat']);
end

result_rk4 = cell(1,2);
result_rk4{1} = load(['/home/hofinger/dev/matlab/workdir/Holger/2018_JR/ConnectomeNoDrivers/2019_06_16_No_DriverWIthMeanY_RK4_double_dt/Connectome1.mat']);
result_rk4{2} = load(['/home/hofinger/dev/matlab/workdir/Holger/2018_JR/ConnectomeNoDrivers/2019_06_16_No_DriverWIthMeanY_RK4/Connectome1.mat']);

result_rk4_less_time = cell(1,2);
for k=1:2
    result_rk4_less_time{k} = load(['/home/hofinger/dev/matlab/workdir/Holger/2018_JR/ConnectomeNoDrivers/2019_06_16_No_DriverWIthMeanY_RK4_less_time/Connectome' num2str(k) '.mat']);
end

%%
result_all = cat(2, result_euler, result_rk4, result_rk4_less_time);

coh_vals = cell(1,length(result_all));
for k=1:length(result_all)
  coh_vals{k} = result_all{k}.simResult.FC{1}(triuIdx);
  
  figure(k+10)
  imagesc(result_all{k}.simResult.FC{1})
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
for k=1:length(result_all)
  plot(result_all{k}.simResult.Y(node_idx,:))
end
hold off;