load('rho_anova_inputs.mat')

tbl = table(groupMatch(:),groupSubjDTI(:),groupSubjEEG(:),metricVal(:),'VariableNames',{'groupMatch','groupSubjDTI','groupSubjEEG','metricVal'});

x = fitlme(tbl,'metricVal~1+groupMatch+(1|groupSubjEEG)+(1|groupSubjDTI)');

std(x.residuals)
ssq(x.residuals)
rssq(x.residuals)^2
rssq(x.residuals)^2