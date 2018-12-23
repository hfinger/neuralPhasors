clear all;

%%
data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['rgast/JansenRit/JR_Paper/' my_foldername]);

results = load(fullfile( path_results, 'all_coh.mat'));

paramSizes = cellfun(@length, results.paramValues);

results.coh = squeeze(reshape(results.all_coh, paramSizes));
results.FC = squeeze(reshape(results.all_FC, paramSizes));
results.corr_SimFC = squeeze(reshape(results.all_corr_SimFC, paramSizes));
results.corr_SimFC = cellfun(@(x) x{1}, results.corr_SimFC);

k1 = cell2mat(results.paramValues{1});
v1 = cell2mat(results.paramValues{2});

%%

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results2 = fullfile(data.resultsdir, ['rgast/JansenRit/JR_Paper/JR_NoDriverOpt2']);

results2 = load(fullfile( path_results2, 'all_coh.mat'));

paramSizes = cellfun(@length, results2.paramValues);

results2.coh = squeeze(reshape(results2.all_coh, paramSizes));
results2.FC = squeeze(reshape(results2.all_FC, paramSizes));
results2.corr_SimFC = squeeze(reshape(results2.all_corr_SimFC, paramSizes));
results2.corr_SimFC = cellfun(@(x) x{1}, results2.corr_SimFC);

k2 = cell2mat(results2.paramValues{1});
v2 = cell2mat(results2.paramValues{2});

%% combine

corr_SimFC_tot = cat(1, results.corr_SimFC, results2.corr_SimFC);
k = cat(2, k1, k2);
v = cat(2, v1, v2);

%%

imagesc(corr_SimFC_tot')
set(gca,'YDir','normal');
set(gca, 'XTick', 1:length(k), 'XTickLabel', k)
set(gca, 'YTick', 1:length(v), 'YTickLabel', v)
xlabel('k')
ylabel('v')
title('correlation to empirical FC')
colorbar;
saveas(gcf,fullfile( 'corr_to_empFC.png'))

%% plot SC and Dist:
p = 1;
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
C = bsxfun(@rdivide,C,sum(C,2));
symSC = (C + C') / 2;

%% evaluate specific parameter set:
empFC = getEmpFC( true, 10, true );
simFC = results2.FC{1,11,1}{1}(3:end,3:end);
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
saveas(gcf,fullfile('C.png'))

figure(3)
imagesc(D)
title('D');
set(gca,'YDir','normal');
colormap(flipud(colormap))
colorbar;
set(gca,'clim',[0, 200])
saveas(gcf,fullfile('D.png'))

figure(4)
imagesc(simFC)
title('simFC');
set(gca,'YDir','normal');
colorbar;
set(gca,'clim',[0, 0.9])
saveas(gcf,fullfile('simFC.png'))

figure(5)
imagesc(empFC)
title('empFC');
set(gca,'YDir','normal');
colorbar;
set(gca,'clim',[0, 0.5])
saveas(gcf,fullfile('empFC.png'))

figure(6)
imagesc(symSC)
title('symSC');
set(gca,'YDir','normal');
colorbar;
set(gca,'clim',[0, 0.6])
saveas(gcf,fullfile('symSC.png'))
