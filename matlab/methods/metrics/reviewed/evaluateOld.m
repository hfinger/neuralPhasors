load('logLearnSC.mat')

%%
logLearnSC = cellfun(@(x) cell2mat(permute(x,[2 3 1])), logLearnSC, 'UniformOutput', false);
logLearnSC = logLearnSC(1:2,:);
finalLearnSC = finalLearnSC(1:2,:);

%%
all1 = cat(3,logLearnSC{:});
all2 = cat(3,finalLearnSC{:});


cell2mat(tmp(1:numRep,:))

% [simFCcov, simFCcor] = sar(logLearnSC{1,1}{end}, 0.65);
% corr( simFCcov(triuIds), empFC(triuIds))

