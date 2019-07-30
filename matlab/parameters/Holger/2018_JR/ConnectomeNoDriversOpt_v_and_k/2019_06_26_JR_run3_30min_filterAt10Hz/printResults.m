data = dataPaths();
[~,my_foldername] = fileparts(pwd);

my_path = ['Holger/2018_JR/ConnectomeNoDriversOpt_v_and_k/' my_foldername];
path_results = fullfile(data.resultsdir, my_path);
path_workdir = fullfile(data.workdir, my_path);

results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

results.freqs = squeeze(reshape(results.all_freqs, paramSizes));
results.coh = squeeze(reshape(results.all_coh, paramSizes));
results.FC = squeeze(reshape(results.all_FC, paramSizes));
results.corr_SimFC = squeeze(reshape(results.all_corr_SimFC, paramSizes));
results.corr_SimFC = cellfun(@(x) x{1}, results.corr_SimFC);

k = cell2mat(results.paramValues{1});
v = cell2mat(results.paramValues{2});

kSelectedIdx = 6;
vSelectedIdx = 6;

kSelected = k(kSelectedIdx);
vSelected = v(vSelectedIdx);

disp(['kSelected = ' num2str(kSelected)]);
disp(['vSelected = ' num2str(vSelected)]);

params.JansenRitConnectomePaper.k = 14; %num2cell([3, 10]); %30; %num2cell(round(22:2:34)); %global connection strength scaling
params.JansenRitConnectomePaper.v = 2.6; %num2cell(2.4:0.1:3.2); %3.2; % velocity [m/s]

disp(['freqSelected = ' num2str(results.freqs{kSelectedIdx,vSelectedIdx})]);

%%
imagesc(results.corr_SimFC')
set(gca, 'XTick', 1:length(k), 'XTickLabel', k)
set(gca, 'YTick', 1:length(v), 'YTickLabel', v)
xlabel('k')
ylabel('v')
set(gca,'YDir','normal')
title('corr with empricial FC')
set(gca,'clim',[0.15; 0.65]);
colorbar
saveas(gcf, fullfile( path_results, 'corr_with_empirical_fc.png'))
saveas(gcf, fullfile( path_results, 'corr_with_empirical_fc.pdf'))

%% plot SC and Dist:
p = 1;
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
C = bsxfun(@rdivide,C,sum(C,2));
symSC = (C + C') / 2;

%% evaluate specific parameter set:
empFC = getEmpFC( true, 10, true );
simFC = results.FC{kSelectedIdx,vSelectedIdx,1}{1}(3:end,3:end);
simFC(logical(eye(33))) = 0;

Idx_mat = triu(ones(size(empFC)),1) == 1;
Idx_mat_offdiag = eye(size(empFC,1)) == 0;

[corr_simFC_empFC,P] = corrcoef(simFC(Idx_mat), empFC(Idx_mat));
disp(['corr_simFC_empFC = ' num2str(corr_simFC_empFC(1,2)) ' p=' num2str(P(1,2))]);

[corr_simFC_D,P] = corrcoef(simFC(Idx_mat), D(Idx_mat));
disp(['corr_simFC_D = ' num2str(corr_simFC_D(1,2)) ' p=' num2str(P(1,2))]);

[corr_empFC_D,P] = corrcoef(empFC(Idx_mat), D(Idx_mat));
disp(['corr_empFC_D = ' num2str(corr_empFC_D(1,2)) ' p=' num2str(P(1,2))]);

[corr_simFC_symSC,P] = corrcoef(simFC(Idx_mat), symSC(Idx_mat));
disp(['corr_simFC_symSC = ' num2str(corr_simFC_symSC(1,2)) ' p=' num2str(P(1,2))]);

[corr_empFC_symSC,P] = corrcoef(empFC(Idx_mat), symSC(Idx_mat));
disp(['corr_empFC_symSC = ' num2str(corr_empFC_symSC(1,2)) ' p=' num2str(P(1,2))]);

[corr_simFC_C,P] = corrcoef(simFC(Idx_mat_offdiag), C(Idx_mat_offdiag));
disp(['corr_simFC_C = ' num2str(corr_simFC_C(1,2)) ' p=' num2str(P(1,2))]);

[corr_empFC_C,P] = corrcoef(empFC(Idx_mat_offdiag), C(Idx_mat_offdiag));
disp(['corr_empFC_C = ' num2str(corr_empFC_C(1,2)) ' p=' num2str(P(1,2))]);

%%
figure(2)
imagesc(C)
title('C');
set(gca,'YDir','normal');
colorbar;
set(gca,'clim',[0, 0.8])
saveas(gcf,fullfile(path_results, 'C.png'))

figure(3)
imagesc(D)
title('D');
set(gca,'YDir','normal');
colormap(flipud(colormap))
colorbar;
set(gca,'clim',[0, 200])
saveas(gcf,fullfile(path_results, 'D.png'))

figure(4)
imagesc(simFC)
title('simFC');
set(gca,'YDir','normal');
colorbar;
set(gca,'clim',[0, 0.9])
saveas(gcf,fullfile(path_results, 'simFC.png'))

figure(5)
imagesc(empFC)
title('empFC');
set(gca,'YDir','normal');
colorbar;
set(gca,'clim',[0, 0.5])
saveas(gcf,fullfile(path_results, 'empFC.png'))

figure(6)
imagesc(symSC)
title('symSC');
set(gca,'YDir','normal');
colorbar;
set(gca,'clim',[0, 0.6])
saveas(gcf,fullfile(path_results, 'symSC.png'))
