clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results1 = fullfile(data.resultsdir, ['Holger/2018_JR/ConnectomeTwoDrivers/AllNodePairs11Hz500uV']);
path_results2 = fullfile(data.resultsdir, ['Holger/2018_JR/ConnectomeTwoDrivers/AllNodePairs11Hz500uVsameSeed']);
results1 = load(fullfile( path_results1, 'all_coh.mat'));
results2 = load(fullfile( path_results2, 'all_coh.mat'));
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ConnectomeTwoDrivers/' my_foldername]);
if ~exist(path_results,'dir')
    mkdir(path_results)
end

results = results1;
for k=1:length(results.all_FC)
    results.all_FC{k}{1} = (results.all_FC{k}{1} + results2.all_FC{k}{1})/2;
end

paramSizes = cellfun(@length, results.paramValues);

all_paths = load(fullfile( data.resultsdir, 'Holger/2018_JR/ConnectomeTwoDrivers', 'paths.mat'));

%% load SC and Dist:
p = 1;
[C,D,F] = getConnectome(1,p,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
C = bsxfun(@rdivide,C,sum(C,2));
symSC = (C + C') / 2;

%%
stim_pair_coh = zeros(length(results.all_FC),1);
for k=1:length(results.all_FC)
    FC = results.all_FC{k}{1}(3:end, 3:end);
    stim_pairs = results.paramComb{2,k};
    stim_pair_coh(k) = FC(stim_pairs(1),stim_pairs(2));
end
stim_pair_coh = squeeze(reshape(stim_pair_coh, paramSizes));
radPerIndex = [0:0.0625:1]*2*pi;

%% we need to remove the last index which is double because 0=2*pi
cohUnique = stim_pair_coh(:,1:16);
radPerIndexUnique = radPerIndex(:,1:16);

[maxVals, maxInd] = max(cohUnique, [], 2);
[minVals, minInd] = min(cohUnique, [], 2);
IPSF = maxVals - minVals;

maxPhaseRad = radPerIndexUnique(maxInd);
minPhaseRad = radPerIndexUnique(minInd);

circHistBinEdges = ([0:0.0625:1] - 0.0625/2)*2*pi;

%% now evaluate shortest paths and PAs:
stimPairs = results.paramValues{2};
shortestPathPerStimPair = zeros(length(stimPairs),1);
distPerStimPair = zeros(length(stimPairs),1);
for k=1:length(stimPairs)
    pair_idx1 = find((all_paths.node_pairs(:,1) == stimPairs{k}(1)) & (all_paths.node_pairs(:,2) == stimPairs{k}(2)));
    pair_idx2 = find((all_paths.node_pairs(:,1) == stimPairs{k}(2)) & (all_paths.node_pairs(:,2) == stimPairs{k}(1)));
    
    shortest1 = all_paths.shortest_path_length(pair_idx1);
    shortest2 = all_paths.shortest_path_length(pair_idx2);
    
    shortestPathPerStimPair(k) = min(shortest1, shortest2);
    
    distPerStimPair(k) = D(stimPairs{k}(1), stimPairs{k}(2));
end

%%
shortestPathUnique = unique(shortestPathPerStimPair);
shortestPathHistCounts = histc(shortestPathPerStimPair, shortestPathUnique);

meanIPSF = zeros(length(shortestPathUnique),1);
stdIPSF = zeros(length(shortestPathUnique),1);
for k=1:length(shortestPathUnique)
    idxs = shortestPathPerStimPair == shortestPathUnique(k);
    tmp = IPSF(idxs);
    meanIPSF(k) = mean(tmp);
    stdIPSF(k) = std(tmp);
end

%% Now evaluate the PA:
paths_grid = cell(33,33);
for k=1:size(all_paths.node_pairs, 1)
    node1 = all_paths.node_pairs(k,1);
    node2 = all_paths.node_pairs(k,2);
    paths_grid{node1, node2} = all_paths.paths_all{k};
end

%%
paths_unique_PAs = cell(length(results.all_FC),1);
paths_unique = cell(length(results.all_FC),1);
paths_twoway = cell(length(results.all_FC),1);
for k=1:length(results.all_FC)
    if mod(k,500)==0
        disp(k)
    end
    FC = results.all_FC{k}{1}(3:end, 3:end);
    stim_pairs = results.paramComb{2,k};
    
    paths1 = paths_grid{stim_pairs(1),stim_pairs(2)};
    paths2 = paths_grid{stim_pairs(2),stim_pairs(1)};
    paths1 = num2cell(paths1,1);
    paths2 = num2cell(paths2,1);
    paths1 = cellfun(@nonzeros, paths1, 'UniformOutput', false);
    paths2 = cellfun(@nonzeros, paths2, 'UniformOutput', false);
    
    % only invert path direction of paths2 to make them compatible and unique:
    paths2 = cellfun(@flipud, paths2, 'UniformOutput', false);
    
    % now make them unique:
    paths_twoway{k} = cat(2, paths1, paths2);
    paths_chars = cellfun(@(x) num2str(x(:)'),paths_twoway{k},'UniformOutput',false);
    [~,unique_idx] = unique(paths_chars);
    paths_unique{k} = paths_twoway{k}(unique_idx);
    
    % now evaluate min(FC) along paths:
    paths_start = cellfun(@(x) x(1:end-1), paths_unique{k}, 'UniformOutput', false);
    paths_final = cellfun(@(x) x(2:end), paths_unique{k}, 'UniformOutput', false);
    
    paths_lin_ind = cellfun(@(x) sub2ind([33, 33], x(1:end-1), x(2:end)), paths_unique{k}, 'UniformOutput', false);
    coh_along_paths = cellfun(@(x) FC(x), paths_lin_ind, 'UniformOutput', false);
    paths_unique_PAs{k} = cellfun(@(x) min(x), coh_along_paths);
end
paths_unique_PAs = squeeze(reshape(paths_unique_PAs, paramSizes));
paths_unique = squeeze(reshape(paths_unique, paramSizes));
paths_twoway = squeeze(reshape(paths_twoway, paramSizes));
paths_unique_lengths = cellfun(@(x) cellfun(@length,x) - 1, paths_unique(:,1), 'UniformOutput', false);
paths_twoway_lengths = cellfun(@(x) cellfun(@length,x) - 1, paths_twoway(:,1), 'UniformOutput', false);

%% filter to include only paths with length <= 5:
filtered_paths_unique_PAs = cell(size(paths_unique_PAs));
filtered_paths_unique = cell(size(paths_unique));
filtered_paths_unique_lengths = cell(size(paths_unique_lengths));
for k=1:length(paths_unique_lengths)
    filtered_idxs = paths_unique_lengths{k} <= 5;
    filtered_paths_unique_idxs{k} = filtered_idxs;
    filtered_paths_unique_PAs{k} = paths_unique_PAs{k}(filtered_idxs);
    filtered_paths_unique{k} = paths_unique{k}(filtered_idxs);
    filtered_paths_unique_lengths{k} = paths_unique_lengths{k}(filtered_idxs);
end

%%
maxPA_idxs = cellfun(@(x) find(x==max(x), 1, 'first'), paths_unique_PAs);
maxPA_idxs = maxPA_idxs(:,1:16);
unique_maxPAs = cellfun(@unique, num2cell(maxPA_idxs,2), 'UniformOutput', false);
num_maxPAs = cellfun(@length, unique_maxPAs);

%%
sorted_PAs = cellfun(@sort, paths_unique_PAs, 'UniformOutput', false);

tmp_max_idx = sub2ind(size(paths_unique_PAs),(1:size(paths_unique_PAs,1))',maxInd);
tmp_min_idx = sub2ind(size(paths_unique_PAs),(1:size(paths_unique_PAs,1))',minInd);

PAsAtMaxSPO = paths_unique_PAs(tmp_max_idx);
PAsAtMinSPO = paths_unique_PAs(tmp_min_idx);

diffPAs = cell(length(PAsAtMaxSPO),1);
for k=1:length(PAsAtMaxSPO)
    diffPAs{k} = PAsAtMaxSPO{k} - PAsAtMinSPO{k};
    [B,sortIdxs] = sort(diffPAs{k}, 2, 'descend');
    diffPAs{k} = diffPAs{k}(sortIdxs);
    PAsAtMaxSPO{k} = PAsAtMaxSPO{k}(sortIdxs);
    PAsAtMinSPO{k} = PAsAtMinSPO{k}(sortIdxs);
end

%% select examplary pairs and plot polarplot with most important PAs for different SPOs:
num_paths = 7;
PAs_per_SPO_maxMean = cell(size(paths_unique_PAs, 1), 1);
PAs_per_SPO_maxMean_node_idxs = cell(size(paths_unique_PAs, 1), 1);
PAs_per_SPO_maxMax = cell(size(paths_unique_PAs, 1), 1);
PAs_per_SPO_maxMaxMin = cell(size(paths_unique_PAs, 1), 1);
for k=1:size(paths_unique_PAs, 1)
    % first select Paths based on their mean PA over all SPOs
    PAs_per_SPO = cat(1, paths_unique_PAs{k,1:16});
    
    % select paths with maxMean:
    meanPA_per_path = mean(PAs_per_SPO,1);
    [~,sortIdxs] = sort(meanPA_per_path, 2, 'descend');
    PAs_per_SPO_maxMean{k} = PAs_per_SPO(:,sortIdxs(1:num_paths));
    PAs_per_SPO_maxMean_node_idxs{k} = paths_unique{k,1}(sortIdxs(1:num_paths));
    
    % select paths with maxMax:
    maxPA_per_path = max(PAs_per_SPO);
    [~,sortIdxs] = sort(maxPA_per_path, 2, 'descend');
    PAs_per_SPO_maxMax{k} = PAs_per_SPO(:,sortIdxs(1:num_paths));
    
    % select paths with maxMaxMin:
    maxminPA_per_path = max(PAs_per_SPO) - min(PAs_per_SPO);
    [~,sortIdxs] = sort(maxminPA_per_path, 2, 'descend');
    PAs_per_SPO_maxMaxMin{k} = PAs_per_SPO(:,sortIdxs(1:num_paths));
end

%%
close all;

%%
figure(1);
hist(IPSF);
title('histogram of IPSF')
xlabel('IPSF')
ylabel('number region pairs')

figure(2);
polarhistogram(maxPhaseRad, circHistBinEdges)
title('polar histogram (counting region pairs) by SPO with max IPSF')

figure(3)
bar(shortestPathUnique, shortestPathHistCounts)
title('number of region pairs per shortest path length')
xlabel('shortest path length')
ylabel('number region pairs')

%% IPSF_vs_shortest_path_length:
[p,tbl,stats] = anova1(IPSF, shortestPathPerStimPair, 'off');
anova_result = ['F = ' num2str(tbl{2, 5}) ', df = ' num2str(tbl{2, 3}) ', p = ' num2str(tbl{2, 6}) ];
print(anova_result)

figure(4); 
clf;
notBoxPlot(IPSF, shortestPathPerStimPair, 'jitter', 0.4)
set(gca,'xlim',[0.5, 6.5])
xlabel('shortest path length')
ylabel('IPSF')
title(anova_result)
saveas(gcf,fullfile(path_results, 'IPSF_vs_shortest_path_length.png'))
saveas(gcf,fullfile(path_results, 'IPSF_vs_shortest_path_length.pdf'))

%% log IPSF_vs_shortest_path_length:
use_shortest_path_lengths = [1, 2, 3, 4, 5, 6];
use_pair_idxs = sum(shortestPathPerStimPair == use_shortest_path_lengths,2) > 0;
[p,tbl,stats] = anova1(log(IPSF(use_pair_idxs)), shortestPathPerStimPair(use_pair_idxs), 'off');
anova_result = ['F(' num2str(tbl{2, 3}) ',' num2str(tbl{3, 3}) ') = ' num2str(tbl{2, 5}) ', p = ' num2str(tbl{2, 6}) ];
disp(['An analysis of variance showed that the effect of shortest path length on log(IPSF) was significant, ' anova_result])

figure(40); 
clf;
notBoxPlot(log(IPSF(use_pair_idxs)), shortestPathPerStimPair(use_pair_idxs), 'jitter', 0.4)
%set(gca,'xlim',[0.5, 6.5])
set(gca,'ylim',[log(0.02) log(1)])
tmp_ticks = [0.02:0.01:0.1 0.2:0.1:1];
tmp_tick_labels = cellfun(@num2str,num2cell(tmp_ticks), 'UniformOutput', false);
tmp_tick_labels([2, 3, 4, 5,6,7,8,10,11,12,13,14,15,16,17]) = {'','','','','','','','','','','','','','',''};
yticks(log(tmp_ticks))
yticklabels(tmp_tick_labels)
xlabel('shortest path length')
ylabel('IPSF')
title(anova_result)
saveas(gcf,fullfile(path_results, 'logIPSF_vs_shortest_path_length.png'))
saveas(gcf,fullfile(path_results, 'logIPSF_vs_shortest_path_length.pdf'))

%% IPSF_vs_fiber_length
figure(5);
clf;
myFit = LinearModel.fit(distPerStimPair,log(IPSF));
h = plot(myFit);
set(h(1),'marker','.')
set(gca,'ylim',[log(0.02) log(1)])
tmp_ticks = [0.02:0.01:0.1 0.2:0.1:1];
tmp_tick_labels = cellfun(@num2str,num2cell(tmp_ticks), 'UniformOutput', false);
tmp_tick_labels([2, 3, 4, 5,6,7,8,10,11,12,13,14,15,16,17]) = {'','','','','','','','','','','','','','',''};
yticks(log(tmp_ticks))
yticklabels(tmp_tick_labels)
[rho,pval] = corr(distPerStimPair, log(IPSF));
title(['corr = ' num2str(rho) ' pval = ' num2str(pval)])
xlabel('average fiber length')
ylabel('IPSF')
saveas(gcf,fullfile(path_results, 'IPSF_vs_fiber_length.png'))
saveas(gcf,fullfile(path_results, 'IPSF_vs_fiber_length.pdf'))

%% IPSF_vs_num_paths no log IPSF
filter_path_lengths = [1; 2; 3; 4; 5];
use_twoway_paths = false;
if use_twoway_paths
    use_path_lengths = paths_twoway_lengths;
else
    use_path_lengths = paths_unique_lengths;
end
num_unique_paths_per_pair = cellfun(@(x) sum(sum(x==filter_path_lengths)), use_path_lengths);
filter_zero_idxs = (num_unique_paths_per_pair ~= 0);
myFit = LinearModel.fit(num_unique_paths_per_pair(filter_zero_idxs),IPSF(filter_zero_idxs));
figure(6);
clf;
h = plot(myFit);
set(h(1),'marker','.')
[rho,pval] = corr(num_unique_paths_per_pair(filter_zero_idxs), IPSF(filter_zero_idxs));
title(['corr = ' num2str(rho) ' pval = ' num2str(pval)])
xlabel('number of paths with length <= 5')
ylabel('IPSF')
saveas(gcf,fullfile(path_results, 'IPSF_vs_num_paths.png'))
saveas(gcf,fullfile(path_results, 'IPSF_vs_num_paths.pdf'))

figure(7); 
clf;
num_unique_paths_per_pair_thresh = num_unique_paths_per_pair(filter_zero_idxs);
num_unique_paths_per_pair_thresh(num_unique_paths_per_pair_thresh > 6) = 6;
notBoxPlot((IPSF(filter_zero_idxs)), num_unique_paths_per_pair_thresh, 'jitter', 0.4)
xlabel('number of paths with length <= 5')
ylabel('IPSF')

%% log IPSF_vs_num_paths
filter_path_lengths = [1; 2; 3; 4; 5];
use_twoway_paths = false;
if use_twoway_paths
    use_path_lengths = paths_twoway_lengths;
else
    use_path_lengths = paths_unique_lengths;
end
num_unique_paths_per_pair = cellfun(@(x) sum(sum(x==filter_path_lengths)), use_path_lengths);
filter_zero_idxs = (num_unique_paths_per_pair ~= 0);
myFit = LinearModel.fit(num_unique_paths_per_pair(filter_zero_idxs),log(IPSF(filter_zero_idxs)));

figure(8);
clf;
h = plot(myFit);
set(h(1),'marker','.')
set(gca,'ylim',[log(0.02) log(1)])
tmp_ticks = [0.02:0.01:0.1 0.2:0.1:1];
tmp_tick_labels = cellfun(@num2str,num2cell(tmp_ticks), 'UniformOutput', false);
tmp_tick_labels([2, 3, 4, 5,6,7,8,10,11,12,13,14,15,16,17]) = {'','','','','','','','','','','','','','',''};
yticks(log(tmp_ticks))
yticklabels(tmp_tick_labels)
[rho,pval] = corr(num_unique_paths_per_pair(filter_zero_idxs), log(IPSF(filter_zero_idxs)));
title(['corr = ' num2str(rho) ' pval = ' num2str(pval)])
xlabel('number of paths with length <= 5')
ylabel('IPSF')
saveas(gcf,fullfile(path_results, 'logIPSF_vs_num_paths.png'))
saveas(gcf,fullfile(path_results, 'logIPSF_vs_num_paths.pdf'))

[p,tbl,stats] = anova1(log(IPSF(filter_zero_idxs)), num_unique_paths_per_pair_thresh, 'off');
anova_result = ['F(' num2str(tbl{2, 3}) ',' num2str(tbl{3, 3}) ') = ' num2str(tbl{2, 5}) ', p = ' num2str(tbl{2, 6}) ];
disp(['An analysis of variance showed that the effect of the number of paths (with path length <= 5, all regions with more than 5 paths were pooled) on log(IPSF) was significant, ' anova_result])

figure(9); 
clf;
num_unique_paths_per_pair_thresh = num_unique_paths_per_pair(filter_zero_idxs);
num_unique_paths_per_pair_thresh(num_unique_paths_per_pair_thresh > 6) = 6;
notBoxPlot(log(IPSF(filter_zero_idxs)), num_unique_paths_per_pair_thresh, 'jitter', 0.4)
yticks(log(tmp_ticks))
yticklabels(tmp_tick_labels)
xlabel('number of paths with length <= 5')
ylabel('IPSF')
title(anova_result)
saveas(gcf,fullfile(path_results, 'logIPSF_vs_num_paths_nobox.png'))
saveas(gcf,fullfile(path_results, 'logIPSF_vs_num_paths_nobox.pdf'))



%% polar plots of max IPSF phases:
figure(10);
clf;
for k=1:6
    subplot(2,3,k)
    idxs = shortestPathPerStimPair == shortestPathUnique(k);
    counts = sum(idxs);
    polarhistogram(maxPhaseRad(idxs), circHistBinEdges)
    thetaticks(0:45:315)
    %rlim([0, 15])
    title(['SP=' num2str(k) ' counts=' num2str(counts)])
end
saveas(gcf,fullfile(path_results, 'polarhistogram_IPSF_phases.png'))
saveas(gcf,fullfile(path_results, 'polarhistogram_IPSF_phases.pdf'))

%% polar plots of max IPSF phases AND filter to include only pairs with IPSF > 0.1:
figure(11);
clf;
for k=1:6
    subplot(3,2,k)
    idxs = ((shortestPathPerStimPair == shortestPathUnique(k)) & (IPSF > 0.1));
    counts = sum(idxs);
    polarhistogram(maxPhaseRad(idxs), circHistBinEdges)
    title(['SP=' num2str(k) ' counts=' num2str(counts)])
end

%%
figure(12);
polarscatter(maxPhaseRad,IPSF)
title('IPSF per SPO');

figure(13);
mean_IPSF_per_polar = zeros(max(maxInd),1);
for k=1:max(maxInd)
    idxs = (maxInd == k);
    mean_IPSF_per_polar(k) = mean(IPSF(idxs));
end
polarplot(radPerIndexUnique,mean_IPSF_per_polar)
title('mean IPSF per SPO');

%% plot PA stats:
figure(14);
clf;
hist(num_maxPAs,min(num_maxPAs):max(num_maxPAs))
xlabel('number of switching pathways')
ylabel('number of region pairs')
title('number of different active paths depending on phase offset')
saveas(gcf,fullfile(path_results, 'histogram_of_PA_counts.png'))

%% plot PAs but filter for real switching without noise switching
figure(15);
clf;
for k=1:6
    subplot(3,2,k)
    idxs = shortestPathPerStimPair == shortestPathUnique(k);
    num_maxPAs_this_sp = num_maxPAs(idxs);
    hist(num_maxPAs_this_sp,min(num_maxPAs):max(num_maxPAs))
    xlabel('number of switching pathways')
    title(['all region pairs with shortest path length = ' num2str(k)])
end

%% use only the pairs with more than 30 pathways:
num_pathways = 20;
use_shortest_path_lengths = [2, 3, 4];
numPathsPerPair = cellfun(@length, diffPAs);
use_pair_idxs = numPathsPerPair >= num_pathways & sum(shortestPathPerStimPair == use_shortest_path_lengths,2) > 0;
diffPAsFiltered = cellfun(@(x) x(1:num_pathways), diffPAs(use_pair_idxs), 'UniformOutput', false);
PAsAtMaxSPOFiltered = cellfun(@(x) x(1:num_pathways), PAsAtMaxSPO(use_pair_idxs), 'UniformOutput', false);
meanDiffPAs = mean(cell2mat(diffPAsFiltered),1);
meanMaxSPOPAs = mean(cell2mat(PAsAtMaxSPOFiltered),1);

figure(16);
clf;
bar(1:num_pathways,meanMaxSPOPAs)
w2 = .8;
hold on;
bar(1:num_pathways,meanDiffPAs,w2,'FaceColor',[1.0 0 0])
hold off;
title('Mean PA at max SPO (blue), diff to min SPO (red)')
xlabel('path index')
ylabel('path activation')
saveas(gcf,fullfile(path_results, 'MeanPAsAtMaxVsMinSPO.png'))

%% use only paths with maximum length 5 and use only the pairs with more than 10 pathways:
use_paths_with_max_lengths = 4;
num_pathways = 10;
use_shortest_path_lengths = [1, 2, 3, 4];
diffPAsFilteredByPathLength = cell(size(diffPAs));
for k=1:length(paths_unique_lengths)
    filtered_idxs = paths_unique_lengths{k} <= use_paths_with_max_lengths;
    diffPAsFilteredByPathLength{k} = diffPAs{k}(filtered_idxs);
end

use_pair_idxs = cellfun(@length, diffPAsFilteredByPathLength) >= num_pathways & sum(shortestPathPerStimPair == use_shortest_path_lengths,2) > 0;
diffPAsFiltered = cellfun(@(x) x(1:num_pathways), diffPAsFilteredByPathLength(use_pair_idxs), 'UniformOutput', false);
PAsAtMaxSPOFiltered = cellfun(@(x) x(1:num_pathways), PAsAtMaxSPO(use_pair_idxs), 'UniformOutput', false);
diffPAsFiltered = cell2mat(diffPAsFiltered);
meanDiffPAs = mean(diffPAsFiltered,1);
meanMaxSPOPAs = mean(cell2mat(PAsAtMaxSPOFiltered),1);

figure(17);
clf;
bar(1:num_pathways,meanMaxSPOPAs)
w2 = .8;
hold on;
bar(1:num_pathways,meanDiffPAs,w2,'FaceColor',[1.0 0 0])
hold off;
title('Mean PA at max SPO (blue), diff to min SPO (red)')
xlabel('path index')
ylabel('path activation')
saveas(gcf,fullfile(path_results, 'MeanPAsAtMaxVsMinSPO_filtered.png'))
saveas(gcf,fullfile(path_results, 'MeanPAsAtMaxVsMinSPO_filtered.pdf'))

figure(18);
h = notBoxPlot(diffPAsFiltered, 'jitter', 0.4);
xlabel('path index (sorted by PA difference)')
ylabel('mean PA difference (per sorting position)')
saveas(gcf,fullfile(path_results, 'MeanPAsAtMaxVsMinSPO_filtered_nobox.png'))
saveas(gcf,fullfile(path_results, 'MeanPAsAtMaxVsMinSPO_filtered_nobox.pdf'))

%% hypothesis test as suggested by Peter:
use_paths_with_max_lengths = 5;
num_pathways = 10;
num_hypo_test = 10000;
use_winzorising = false;
use_shortest_path_lengths = 1:use_paths_with_max_lengths;
diffPAsFilteredByPathLength = cell(size(diffPAs));
for k=1:length(paths_unique_lengths)
    filtered_idxs = paths_unique_lengths{k} <= use_paths_with_max_lengths;
    diffPAsFilteredByPathLength{k} = diffPAs{k}(filtered_idxs);
end
pValues = cell(size(diffPAsFilteredByPathLength));
for k=1:length(paths_unique_lengths)
    tmp = diffPAsFilteredByPathLength{k};
    if length(tmp) >= num_pathways
        if use_winzorising
            meanTmp = mean(tmp(2:end-1));
            stdTmp = std(tmp(2:end-1));
        else
            meanTmp = mean(tmp);
            stdTmp = std(tmp);
        end
        testingDist = meanTmp + stdTmp * randn(length(tmp), num_hypo_test);
        maxOfTestingDist = max(testingDist,[],1);
        pValues{k} = sum(logical(sum(maxOfTestingDist > tmp(1), 1))) / num_hypo_test;
    end
end

use_pair_idxs = cellfun(@length, diffPAsFilteredByPathLength) >= num_pathways & sum(shortestPathPerStimPair == use_shortest_path_lengths,2) > 0;
diffPAsFilteredByPathLengthFilteredPairs = diffPAsFilteredByPathLength(use_pair_idxs);
pValues = pValues(use_pair_idxs);
pValues = cell2mat(pValues);
pValues(pValues < (1/num_hypo_test)) = 1 / num_hypo_test;
chi_vals = -2.*log(pValues);
fisher_group_pval = 1 - chi2cdf(sum(chi_vals),2*length(pValues));
geom_pval = geomean(pValues);

figure(19);
hist(pValues, 40)
title(['geomean=' num2str(geom_pval) ' fisher p=' num2str(fisher_group_pval)])
xlabel('p')
ylabel('number of region pairs')
saveas(gcf,fullfile(path_results, 'hist_of_pvalues_of_PAdiff.png'))

%% now the same as above but with normalizing:
use_paths_with_min_lengths = 5;
num_pathways = 10;
use_shortest_path_lengths = [2, 3, 4, 5];

diffPAsFilteredNormed = cell(length(PAsAtMaxSPO),1);
PAsAtMinSPOFilteredNormed = cell(length(PAsAtMaxSPO),1);
PAsAtMaxSPOFilteredNormed = cell(length(PAsAtMaxSPO),1);
for k=1:length(PAsAtMaxSPO)
    
    filtered_idxs = paths_unique_lengths{k} <= use_paths_with_min_lengths;
    
    tmpPAMaxSPO = PAsAtMaxSPO{k}(filtered_idxs);
    tmpPAMinSPO = PAsAtMinSPO{k}(filtered_idxs);
    
    tmpPAMaxSPO = tmpPAMaxSPO / sum(tmpPAMaxSPO);
    tmpPAMinSPO = tmpPAMinSPO / sum(tmpPAMinSPO);
    
    diffPAsnormed = tmpPAMaxSPO - tmpPAMinSPO;
    [B,sortIdxs] = sort(diffPAsnormed, 2, 'descend');
    diffPAsFilteredNormed{k} = diffPAsnormed(sortIdxs);
    PAsAtMaxSPOFilteredNormed{k} = tmpPAMaxSPO(sortIdxs);
    PAsAtMinSPOFilteredNormed{k} = tmpPAMinSPO(sortIdxs);
end

use_pair_idxs = cellfun(@length, diffPAsFilteredNormed) >= num_pathways & sum(shortestPathPerStimPair == use_shortest_path_lengths,2) > 0;
diffPAsFilteredNormed = cellfun(@(x) x(1:num_pathways), diffPAsFilteredNormed(use_pair_idxs), 'UniformOutput', false);
PAsAtMaxSPOFiltered = cellfun(@(x) x(1:num_pathways), PAsAtMaxSPO(use_pair_idxs), 'UniformOutput', false);
diffPAsFilteredNormed = cell2mat(diffPAsFilteredNormed);
meanDiffPAsFilteredNormed = mean(diffPAsFilteredNormed,1);
meanMaxSPOPAsFilteredNormed = mean(cell2mat(PAsAtMaxSPOFiltered),1);

figure(20);
clf;
bar(1:num_pathways,meanMaxSPOPAsFilteredNormed)
w2 = .8;
hold on;
bar(1:num_pathways,meanDiffPAsFilteredNormed,w2,'FaceColor',[1.0 0 0])
hold off;
title('Mean PA at max SPO (blue), diff to min SPO (red)')
xlabel('path index')
ylabel('path activation')


%% compare most active path with second most active path that has no overlap with first:
use_max_path_length = 5;
PAs_of_collected_paths_cell = cell(1,1,length(paths_unique_lengths));
first_path_idxs_cell = cell(1,1,length(paths_unique_lengths));
first_path_PAs_cell = cell(1,1,length(paths_unique_lengths));
second_path_idxs_cell = cell(1,1,length(paths_unique_lengths));
PAs_of_collected_third_paths = cell(1,1,length(paths_unique_lengths));
IPSFs_of_collected_paths = cell(1,1,length(paths_unique_lengths));
maxSpoIdx_of_collected_paths = cell(1,1,length(paths_unique_lengths));
cohUnique_of_collected_paths = cell(1,1,length(paths_unique_lengths));
for k=1:length(paths_unique_lengths)
    use_paths_idx = (use_max_path_length >= paths_unique_lengths{k});
    this_PAs = cell2mat(paths_unique_PAs(k,:)');
    avg_PAs = mean(this_PAs,1);
    [~, sidx] = sort(avg_PAs, 2, 'descend');
    
    first_path_idx = sidx(1);
    first_path_idxs_cell{1,1,k} = first_path_idx;
    first_path_PA = this_PAs(:, first_path_idx);
    nodes_on_first_path = paths_unique{k,1}{first_path_idx}(2:end-1);
    first_path_PAs_cell{1,1,k} = first_path_PA;
    
    % now search for second best without overlap:
    for j=2:length(sidx)
        second_path_idx = sidx(j);
        nodes_on_second_path = paths_unique{k,1}{second_path_idx}(2:end-1);
        if isempty(intersect(nodes_on_first_path, nodes_on_second_path))
            % found disjoint path:
            second_path_PA = this_PAs(:, second_path_idx);
            PAs_of_collected_paths_cell{1,1,k} = [first_path_PA, second_path_PA];
            second_path_idxs_cell{1,1,k} = second_path_idx;
            IPSFs_of_collected_paths{1,1,k} = IPSF(k);
            maxSpoIdx_of_collected_paths{1,1,k} = maxInd(k);
            cohUnique_of_collected_paths{1,1,k} = cohUnique(k,:);
            
            % now search for third best without overlap:
            for l=(j+1):length(sidx)
                third_path_idx = sidx(l);
                nodes_on_third_path = paths_unique{k,1}{third_path_idx}(2:end-1);
                if isempty(intersect(nodes_on_first_path, nodes_on_third_path)) && isempty(intersect(nodes_on_second_path, nodes_on_third_path))
                    % found third disjoint path:
                    third_path_PA = this_PAs(:, third_path_idx);
                    PAs_of_collected_third_paths{1,1,k} = [first_path_PA, third_path_PA];
                    break;
                end
            end
            
            break;
        end
    end
end

PAs_of_collected_paths = cell2mat(PAs_of_collected_paths_cell);
PAs_of_collected_third_paths = cell2mat(PAs_of_collected_third_paths);
IPSFs_of_collected_paths = cell2mat(IPSFs_of_collected_paths);
maxSpoIdx_of_collected_paths = cell2mat(maxSpoIdx_of_collected_paths);
cohUnique_of_collected_paths = squeeze(cell2mat(cohUnique_of_collected_paths));
cohUnique_of_collected_paths_ztrafo = bsxfun(@minus, cohUnique_of_collected_paths, mean(cohUnique_of_collected_paths,1));
cohUnique_of_collected_paths_ztrafo = bsxfun(@rdivide, cohUnique_of_collected_paths_ztrafo, std(cohUnique_of_collected_paths_ztrafo,1));

maxSpoRad_of_collected_paths = radPerIndex(maxSpoIdx_of_collected_paths);

unit_circle = exp(1i * radPerIndex(1:16))';
circMean = squeeze(mean( bsxfun(@times, PAs_of_collected_paths(1:16, :, :), unit_circle), 1));
angleWithMostPA = mod( angle(circMean) ,2*pi);
absWithMostPA = abs(circMean);

[~, maxIndPAs1] = max(PAs_of_collected_paths(1:16, 1, :), [], 1);
[~, maxIndPAs2] = max(PAs_of_collected_paths(1:16, 2, :), [], 1);
maxSPO_PAs1 = radPerIndex(maxIndPAs1);
maxSPO_PAs2 = radPerIndex(maxIndPAs2);
angleWithMostPA_using_max = cat(1, maxSPO_PAs1, maxSPO_PAs2);

cohZ_at_maxIndPAs1 = diag(cohUnique_of_collected_paths_ztrafo(maxIndPAs1,:));
cohZ_at_maxIndPAs2 = diag(cohUnique_of_collected_paths_ztrafo(maxIndPAs2,:));

angleWithMostPA_relToMaxSpo = mod( angleWithMostPA_using_max - maxSpoRad_of_collected_paths, 2*pi);

angleWithMostPA_relToMaxSpo_histCounts1 = histcounts(angleWithMostPA_relToMaxSpo(1,:), circHistBinEdges);
angleWithMostPA_relToMaxSpo_histCounts2 = histcounts(angleWithMostPA_relToMaxSpo(2,:), circHistBinEdges);

expected = sum(angleWithMostPA_relToMaxSpo_histCounts1) / 16;
chi_square1 = sum((expected - angleWithMostPA_relToMaxSpo_histCounts1).^2 / expected);
chi_square2 = sum((expected - angleWithMostPA_relToMaxSpo_histCounts2).^2 / expected);
p1 = chi2cdf(chi_square1,16-1,'upper');
p2 = chi2cdf(chi_square2,16-1,'upper');

figure(21);
clf;
polarhistogram(angleWithMostPA_relToMaxSpo(1,:),circHistBinEdges)
pax = gca; pax.ThetaAxisUnits = 'radians';
thetaticks([0, pi/2, pi, 3*pi/2])
rlim([0, 120])
rticks([0, 120])
saveas(gcf,fullfile(path_results, ['angleWithMostPA_relToMaxSpo_Path1.png']))
saveas(gcf,fullfile(path_results, ['angleWithMostPA_relToMaxSpo_Path1.pdf']))

figure(22);
clf;
polarhistogram(angleWithMostPA_relToMaxSpo(2,:),circHistBinEdges)
pax = gca; pax.ThetaAxisUnits = 'radians';
thetaticks([0, pi/2, pi, 3*pi/2])
rlim([0, 120])
rticks([0, 120])
saveas(gcf,fullfile(path_results, ['angleWithMostPA_relToMaxSpo_Path2.png']))
saveas(gcf,fullfile(path_results, ['angleWithMostPA_relToMaxSpo_Path2.pdf']))

diffAngle = mod( angleWithMostPA_relToMaxSpo(1,:) - angleWithMostPA_relToMaxSpo(2,:), 2*pi);

where_PA1_larger_than_PA2 = squeeze(PAs_of_collected_paths(1:16, 1, :) > PAs_of_collected_paths(1:16, 2, :));
fraction_of_switchers = mean(sum(where_PA1_larger_than_PA2, 1)<14);

figure(77)
hist(mean(where_PA1_larger_than_PA2, 1));
saveas(gcf,fullfile(path_results, ['hist_of_fraction_SPO_where_PA1_larger_than_PA2.pdf']))

[switching_index] = pathway_switching_index(PAs_of_collected_paths(1:16, 1, :),PAs_of_collected_paths(1:16, 2, :));

hist_edges = exp(-5:0.5:-1);
hist_edges = [-hist_edges(end:-1:1), 0, hist_edges];
hist_pair_counts = histcounts(switching_index, hist_edges);

figure(23);
clf;
bar(hist_pair_counts);
set(gca,'xtick',[0.5:length(hist_pair_counts)]);
set(gca,'xticklabels',round(hist_edges*1000)/1000)
xlabel('pathway switching index (PSI)')
ylabel('number node pairs')

%%
figure(24);
clf;
histogram(switching_index, -0.5:0.02:0.5)
xlim([-0.5 0.5])
set(gca,'xtick',[-0.5, 0, 0.5])
set(gca,'ytick',[0, 40])
xlabel('pathway switching index (PSI)')
ylabel('number node pairs')
hold on;
plot([0, 0], [0, 40], 'LineWidth', 2)
hold off;
saveas(gcf,fullfile(path_results, ['hist_PSI.pdf']))

[h,p,ci,stats] = ttest(switching_index);

first_path_PAs = squeeze(cell2mat(first_path_PAs_cell));
first_path_PAs = first_path_PAs(1:16, :);
all_first_path_SPO_selectivity = max(first_path_PAs, [], 1) - min(first_path_PAs, [], 1);

figure(88);
clf;
edges = 0:0.05:0.8;
edges(end) = 0.81;
histogram(all_first_path_SPO_selectivity, edges)
xlim([0 0.8])
set(gca,'xtick',[0 0.4 0.8])
set(gca,'ytick',[0, 180])
xlabel('pathway phase selectivity (PPS)')
ylabel('number node pairs')
saveas(gcf,fullfile(path_results, ['hist_PPS.pdf']))

%%
showPSI = true;
plotOnlyNonOverlapping = true;

figure(25);
set(gcf, 'Position',  [100, 100, 1600, 800])
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
num_per_grpup = 8;
plotDataX = radPerIndexUnique';
plotDataX = [plotDataX; plotDataX(1,:)];
for rng_seed = 1:8
    clf;
    plot_counter = 1;
    rng(rng_seed);
    for plot_pairs_with_shortest_path_length = 1:3
        idxs = find(shortestPathPerStimPair == plot_pairs_with_shortest_path_length);
        idxs = idxs(randperm(length(idxs),num_per_grpup));
        for repeat=1:num_per_grpup
            global_pair_idx = idxs(repeat);
            subplot(3,num_per_grpup,plot_counter)
            plot_counter = plot_counter + 1;
            
            if plotOnlyNonOverlapping
                
                this_PAs = cell2mat(paths_unique_PAs(global_pair_idx,:)');
                this_PAs = [this_PAs(1:16, :); this_PAs(1, :)];
                
                remaining_path_idxs = setdiff(1:size(this_PAs,2), [first_path_idx, second_path_idx]);
                rest_PAs = this_PAs(:,remaining_path_idxs);
                [~, rest_PAs_sort_idxs] = sort(mean(rest_PAs(1:16,:),1),'descend');
                
                
                % first plot remaining and only at last the strongest
                polarplot(plotDataX, rest_PAs(:,rest_PAs_sort_idxs(1)), 'Color', [0.4660, 0.6740, 0.1880]	)
                hold on;
                polarplot(plotDataX, rest_PAs(:,rest_PAs_sort_idxs(2)), 'Color', [0.4660, 0.6740, 0.1880]	)
                polarplot(plotDataX, rest_PAs(:,rest_PAs_sort_idxs(3)), 'Color', [0.4660, 0.6740, 0.1880]	)
                
                second_path_idx = second_path_idxs_cell{1,1,global_pair_idx};
                second_path_PAs = this_PAs(:, second_path_idx);
                if size(second_path_PAs,2)>0
                    polarplot(plotDataX, second_path_PAs, 'Color', [0.6350, 0.0780, 0.1840])
                end
                
                first_path_idx = first_path_idxs_cell{1,1,global_pair_idx};
                first_path_PAs = this_PAs(:, first_path_idx);
                polarplot(plotDataX, first_path_PAs, 'Color', [0, 0.4470, 0.7410])
                
                hold off;
            else
                plotDataY = PAs_per_SPO_maxMean{global_pair_idx};
                %plotDataY = PAs_per_SPO_maxMax{idxs(repeat)};
                %plotDataY = PAs_per_SPO_maxMaxMin{idxs(repeat)};
                
                plotDataY = [plotDataY; plotDataY(1,:)];
                polarplot(plotDataX, plotDataY)
            end
            
            pax = gca; pax.ThetaAxisUnits = 'radians';
            %thetaticks(0:90:315)
            thetaticks([0, pi/2, pi, 3*pi/2])
            
            hold on;
            
            %plotDataYcoh = [cohUnique(idxs(repeat),:), cohUnique(idxs(repeat),1)]';
            %polarplot(plotDataX, plotDataYcoh, 'color', 'k', 'LineStyle', ':')
            
            tmp_rlim = rlim();
            tmp_rlim(1) = 0;
            rlim(tmp_rlim)
            rticks(tmp_rlim(2))
            
            %polarplot([maxPhaseRad(idxs(repeat)), maxPhaseRad(idxs(repeat))], [0, 1])
            
            hold off;
            
            if rng_seed==3 && plot_pairs_with_shortest_path_length==1 && repeat==1
                surfacePlotIdx{1} = idxs(repeat);
            end
            if rng_seed==5 && plot_pairs_with_shortest_path_length==2 && repeat==6
                surfacePlotIdx{2} = idxs(repeat);
            end
            if rng_seed==8 && plot_pairs_with_shortest_path_length==2 && repeat==1
                surfacePlotIdx{3} = idxs(repeat);
            end
            
            if showPSI
                if size(PAs_of_collected_paths_cell{idxs(repeat)},1)>0
                    switching_index = pathway_switching_index(PAs_of_collected_paths_cell{idxs(repeat)}(1:16, 1),PAs_of_collected_paths_cell{idxs(repeat)}(1:16, 2));
                    PPS = all_first_path_SPO_selectivity(idxs(repeat));
                    title(['PSI=' num2str(round(switching_index*100)/100) ' PPS=' num2str(round(PPS*100)/100)])
                end
            else
                if repeat==round(num_per_grpup/2)
                    title(['PA per SPO for 5 randomly selected pairs with shortest path length ' num2str(plot_pairs_with_shortest_path_length) '. Colors correspond to different paths connecting the pair.'])
                end
            end
            
        end
    end
    saveas(gcf,fullfile(path_results, ['PA_per_SPO_for_random_pairs_rng' num2str(rng_seed) '.pdf']))
end

%% surface plots:
[C,D,F] = getConnectome(1,1,0.1,1); % get Connectivity and Delay matrix of Connectome
threshold = 0.1;
C(C < threshold) = 0;
C = bsxfun(@rdivide,C,sum(C,2));
C = C + C';
edgeVals = [C zeros(33,33); zeros(33,33) zeros(33,33)];

all_FCs = cellfun(@(x) x{1}, reshape(results.all_FC, paramSizes), 'UniformOutput', false);

for selected_pair_idx=3 %1:4
    use_pair_idx = surfacePlotIdx{selected_pair_idx};
    my_FC = cell2mat(all_FCs(1,use_pair_idx,:));

    my_FC2 = my_FC(3:end, 3:end, 1:16);

    node_coh_with_driver = squeeze(my_FC(3:end, 1, 1:16));

    my_FC2_whereConn = zeros(size(my_FC2));
    for k=1:16
        tmp = my_FC2(:,:,k);
        tmp(C==0) = 0;
        my_FC2_whereConn(:,:,k) = tmp;
    end

    target_stim_pairs = results.paramComb{2,use_pair_idx};
    coh_between_stimulated_pairs = squeeze(my_FC2(target_stim_pairs(1),target_stim_pairs(2), :));

    for spo_idx = 1:16 %[3, 13]

        % now create a FC conn mat with only connections on given paths:
        my_FC3_wherePath = zeros(size(my_FC2));
        tmp = PAs_per_SPO_maxMean_node_idxs{use_pair_idx};
        tmp_PA = PAs_per_SPO_maxMean{use_pair_idx};
        [~, sort_idx] = sort(tmp_PA(spo_idx,:), 'descend');
        for k=1:1%length(surfacePlotIdx3_node_idxs)
            nodes_on_path = tmp{sort_idx(k)};
            for m=1:length(nodes_on_path)-1
                start_node = nodes_on_path(m);
                end_node = nodes_on_path(m+1);
                my_FC3_wherePath(start_node,end_node,:) = my_FC2_whereConn(start_node,end_node,:);
            end
        end

        edgeWidth = zeros(33,33);
        edgeWidth(my_FC2_whereConn(:,:,spo_idx)>0) = 0.3;
        edgeWidth(my_FC3_wherePath(:,:,spo_idx)>0) = 0.9;
        edgeWidth = [edgeWidth(:,:) zeros(33,33); zeros(33,33) zeros(33,33)];
        edgeColor = [my_FC2_whereConn(:,:,spo_idx) zeros(33,33); zeros(33,33) zeros(33,33)];

        figure(25);
        nodeVals = 0.0001 * zeros(66,1);
        nodeVals(target_stim_pairs) = 0.5;
        nodeVals(1:33) = node_coh_with_driver(:, spo_idx);

        nodeVals(1:33) = 0.2;
        nodeVals(target_stim_pairs) = 0.5;

        nodeScale = 40;
        nodeMinSize = 0;
        nodeRange = [min(min(min(nodeVals))),max(max(max(nodeVals)))];
        edgeMin = 0;
        edgeAlphas = 1;
        surfaceVisibility = 0.1;
        rescaleEdgeVals = 2;
        cmap = colormap(flipud(hot));
        plotBrainConnectivityNew(nodeVals,edgeWidth,nodeMinSize,nodeScale,nodeRange,edgeColor,edgeMin,edgeAlphas,surfaceVisibility,rescaleEdgeVals,1,cmap);
        view(260,15);
        colorbar;
        pause(0.1);

        saveas(gcf,fullfile(path_results, ['connectome_pathways_spo' num2str(spo_idx) '.png']))
        saveas(gcf,fullfile(path_results, ['connectome_pathways_spo' num2str(spo_idx) '.pdf']))
    end
end