clear all;

data = dataPaths();
[~,my_foldername] = fileparts(pwd);
path_results = fullfile(data.resultsdir, ['Holger/2018_JR/ConnectomeTwoDrivers/' my_foldername]);
results = load(fullfile( path_results, 'all_coh.mat'));

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
PAs = cell(length(results.all_FC),1);
PAs_unique = cell(length(results.all_FC),1);
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
    paths_this_sim = cat(2, paths1, paths2);
    paths_chars = cellfun(@(x) num2str(x(:)'),paths_this_sim,'UniformOutput',false);
    [~,unique_idx] = unique(paths_chars);
    PAs_unique{k} = paths_this_sim(unique_idx);
    
    % now evaluate min(FC) along paths:
    paths_start = cellfun(@(x) x(1:end-1), PAs_unique{k}, 'UniformOutput', false);
    paths_final = cellfun(@(x) x(2:end), PAs_unique{k}, 'UniformOutput', false);
    
    paths_lin_ind = cellfun(@(x) sub2ind([33, 33], x(1:end-1), x(2:end)), PAs_unique{k}, 'UniformOutput', false);
    coh_along_paths = cellfun(@(x) FC(x), paths_lin_ind, 'UniformOutput', false);
    PAs{k} = cellfun(@(x) min(x), coh_along_paths);
end
PAs = squeeze(reshape(PAs, paramSizes));
PAs_unique = squeeze(reshape(PAs_unique, paramSizes));

%%
maxPA_idxs = cellfun(@(x) find(x==max(x), 1, 'first'), PAs);
maxPA_idxs = maxPA_idxs(:,1:16);
unique_maxPAs = cellfun(@unique, num2cell(maxPA_idxs,2), 'UniformOutput', false);
num_maxPAs = cellfun(@length, unique_maxPAs);

%%
sorted_PAs = cellfun(@sort, PAs, 'UniformOutput', false);

tmp_max_idx = sub2ind(size(PAs),(1:size(PAs,1))',maxInd);
tmp_min_idx = sub2ind(size(PAs),(1:size(PAs,1))',minInd);

PAsAtMaxSPO = PAs(tmp_max_idx);
PAsAtMinSPO = PAs(tmp_min_idx);

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
PAs_per_SPO_maxMean = cell(size(PAs, 1), 1);
PAs_per_SPO_maxMax = cell(size(PAs, 1), 1);
PAs_per_SPO_maxMaxMin = cell(size(PAs, 1), 1);
for k=1:size(PAs, 1)
    % first select Paths based on their mean PA over all SPOs
    PAs_per_SPO = cat(1, PAs{k,1:16});
    
    % select paths with maxMean:
    meanPA_per_path = mean(PAs_per_SPO,1);
    [~,sortIdxs] = sort(meanPA_per_path, 2, 'descend');
    PAs_per_SPO_maxMean{k} = PAs_per_SPO(:,sortIdxs(1:num_paths));
    
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

%% not box blot:
figure(4); 
clf;
notBoxPlot(IPSF, shortestPathPerStimPair, 'jitter', 0.4)
set(gca,'xlim',[0.5, 6.5])
xlabel('shortest path length')
ylabel('IPSF')
title('IPSF vs shortest path length')
saveas(gcf,fullfile('IPSF_vs_shortest_path_length.png'))
saveas(gcf,fullfile( path_results, 'IPSF_vs_shortest_path_length.png'))
saveas(gcf,fullfile( 'IPSF_vs_shortest_path_length.pdf'))

%%
figure(16);
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
saveas(gcf,fullfile('IPSF_vs_fiber_length.png'))
saveas(gcf,fullfile( 'IPSF_vs_fiber_length.pdf'))

%%
figure(17);
clf;
filter_path_lengths = [1; 2; 3; 4];
path_lengths = cellfun(@(x) cellfun(@length,x), PAs_unique(:,1), 'UniformOutput', false);
num_unique_paths_per_pair = cellfun(@(x) sum(sum(x==filter_path_lengths)), path_lengths);
myFit = LinearModel.fit(num_unique_paths_per_pair,log(IPSF));
h = plot(myFit);
set(h(1),'marker','.')
set(gca,'ylim',[log(0.02) log(1)])
tmp_ticks = [0.02:0.01:0.1 0.2:0.1:1];
tmp_tick_labels = cellfun(@num2str,num2cell(tmp_ticks), 'UniformOutput', false);
tmp_tick_labels([2, 3, 4, 5,6,7,8,10,11,12,13,14,15,16,17]) = {'','','','','','','','','','','','','','',''};
yticks(log(tmp_ticks))
yticklabels(tmp_tick_labels)
[rho,pval] = corr(num_unique_paths_per_pair, log(IPSF));
title(['corr = ' num2str(rho) ' pval = ' num2str(pval)])
xlabel('number of paths with length <= 4')
ylabel('IPSF')
saveas(gcf,fullfile('IPSF_vs_num_paths.png'))
saveas(gcf,fullfile( 'IPSF_vs_num_paths.pdf'))

figure(18); 
clf;
notBoxPlot(log(IPSF), num_unique_paths_per_pair, 'jitter', 0.4)
xlabel('number of paths with length <= 4')
ylabel('IPSF')

%% polar plots of max IPSF phases:
figure(5);
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
saveas(gcf,fullfile('polarhistogram_IPSF_phases.png'))
saveas(gcf,fullfile('polarhistogram_IPSF_phases.pdf'))

%% polar plots of max IPSF phases AND filter to include only pairs with IPSF > 0.1:
figure(6);
clf;
for k=1:6
    subplot(3,2,k)
    idxs = ((shortestPathPerStimPair == shortestPathUnique(k)) & (IPSF > 0.1));
    counts = sum(idxs);
    polarhistogram(maxPhaseRad(idxs), circHistBinEdges)
    title(['SP=' num2str(k) ' counts=' num2str(counts)])
end

%%
figure(7);
polarscatter(maxPhaseRad,IPSF)
title('IPSF per SPO');

figure(8);
mean_IPSF_per_polar = zeros(max(maxInd),1);
for k=1:max(maxInd)
    idxs = (maxInd == k);
    mean_IPSF_per_polar(k) = mean(IPSF(idxs));
end
polarplot(radPerIndexUnique,mean_IPSF_per_polar)
title('mean IPSF per SPO');

%% plot PA stats:
figure(9);
clf;
hist(num_maxPAs,min(num_maxPAs):max(num_maxPAs))
xlabel('number of switching pathways')
ylabel('number of region pairs')
title('number of different active paths depending on phase offset')
saveas(gcf,fullfile('histogram_of_PA_counts.png'))

%% plot PAs but filter for real switching without noise switching
figure(10);
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

figure(11);
clf;
bar(1:num_pathways,meanMaxSPOPAs)
w2 = .8;
hold on;
bar(1:num_pathways,meanDiffPAs,w2,'FaceColor',[1.0 0 0])
hold off;
title('Mean PA at max SPO (blue), diff to min SPO (red)')
xlabel('path index')
ylabel('path activation')
saveas(gcf,fullfile('MeanPAsAtMaxVsMinSPO.png'))

%%
figure(12);
clf;
plot_counter = 1;
for plot_pairs_with_shortest_path_length = 1:3
    idxs = find(shortestPathPerStimPair == plot_pairs_with_shortest_path_length);
    idxs = idxs(randperm(length(idxs),5));
    for repeat=1:5
        
        subplot(3,5,plot_counter)
        plotDataX = radPerIndexUnique';
        %plotDataY = PAs_per_SPO_maxMean{idxs(repeat)};
        %plotDataY = PAs_per_SPO_maxMax{idxs(repeat)};
        plotDataY = PAs_per_SPO_maxMaxMin{idxs(repeat)};
        plotDataX = [plotDataX; plotDataX(1,:)];
        plotDataY = [plotDataY; plotDataY(1,:)];
        polarplot(plotDataX, plotDataY)
        thetaticks(0:45:315)
        rlim([0, 0.9])
        if repeat==3
            title(['PA per SPO for 5 randomly selected pairs with shortest path length ' num2str(plot_pairs_with_shortest_path_length) '. Colors correspond to different paths connecting the pair.'])
        end
        plot_counter = plot_counter + 1;
    end
end
saveas(gcf,fullfile('PA_per_SPO_for_random_pairs.png'))