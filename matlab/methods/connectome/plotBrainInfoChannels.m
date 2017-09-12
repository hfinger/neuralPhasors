function [ results ] = plotBrainInfoChannels( paths, dataStruct, k, dv )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(dataStruct);
n = length(fnames);
paths_ordered = cell(n,k);
pathlengths_ordered = zeros(n,k);
indices = zeros(n,k);
poColl = zeros(n,1);
dvColl = zeros(n,1);

for f=1:length(fnames)
    
    data = dataStruct.(fnames{f});
    targets = data.drivPos;
    targets = targets - length(targets);
    Coh = 1 ./ data.Coherence{1,1}(length(targets)+1:end,length(targets)+1:end);
    poColl(f) = data.drivPO;
    if ~isempty(dv)
        dvColl(f) = data.(dv);
    end
    
    % compute pathlengths of SWPs based on Coherence matrix
    pathlengths = zeros(1,k);
    for i=1:k
        for j=1:length(paths{1,i})-1
            pathlengths(i) = pathlengths(i) + Coh(paths{1,i}(j),paths{1,i}(j+1));
        end
    end

    % order the pathlengths and get indices
    [pathlengths_ordered(f,:),idx] = sort(pathlengths);

    % order paths according to indices
    for i=1:length(idx)
        paths_ordered{f,i} = paths{1,idx(i)};
    end
    indices(f,:) = idx;

end

dvals = unique(dvColl);
for d=1:length(dvals)
    
    idx1 = dvColl == dvals(d);
    dvColl_tmp = dvColl(idx1);
    poColl_tmp = poColl(idx1);
    paths_ordered_tmp = paths_ordered(idx1,:);
    indices_tmp = indices(idx1,:);
    fnames_tmp = fnames(idx1);
    
    [poColl_tmp, idx2] = sort(poColl_tmp);
    dvColl_tmp = dvColl_tmp(idx2);
    paths_ordered_tmp = paths_ordered_tmp(idx2,:);
    indices_tmp = indices_tmp(idx2,:);
    fnames_tmp = fnames_tmp(idx2);
    
    indVals = unique(indices_tmp(:,1));
    edgeVals = zeros(size(Coh,1),size(Coh,2),length(indVals));
    for v=1:length(indVals)
        
        idx3 = indices_tmp(:,1) == indVals(v);
        fnames_tmp2 = fnames_tmp(idx3);
        paths_ordered_tmp2 = paths_ordered_tmp(idx3);
        phaseOffsets{v} = poColl_tmp(idx3)';
        
        % get coherence of simulation to plot
        Coh = dataStruct.(fnames_tmp2{1}).Coherence{1,1}(length(targets)+1:end,length(targets)+1:end);
        
        % get node and edge values for path plotting
        
        for i=1:1
            for j=1:length(paths_ordered_tmp2{1,i})-1
                edgeVals(paths_ordered_tmp2{1,i}(j),paths_ordered_tmp2{1,i}(j+1),v) = 1;
            end
        end
        
    end
    
    nodeVals = zeros(1,size(Coh,1));
    nodeVals(targets) = 1;
    edgeVals = edgeVals * 10;

    % plot nodes and edges in brain
    addpath('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/wetransfer-b16a3e')
    edgeMin = 0.0001; % cut-off value for edgeVals (values below won't be plotted)
    plotBrainConnectivity(nodeVals,edgeVals,edgeMin,1,0)
    if length(dvals) > 1
        set(gcf,'name',['Information Channels in the brain for different Phase Offsets with ',dv,' = ',num2str(dvals(d))])
    end
    
    % include legend for information channels
    lines = findobj(gca,'Type','line');
    legendStr = cell(1,length(phaseOffsets));
    poDiff_standard = poColl_tmp(2) - poColl_tmp(1);
    for i=1:length(phaseOffsets)
        if length(phaseOffsets{1,i}) > 2
            poDiffs = [0,diff(phaseOffsets{1,i})];
            jumps = [1,find(poDiffs > 1.5*poDiff_standard),length(poDiffs)+1];
            poStr = [];
            for j=1:length(jumps)-1
                poStr = [poStr, num2str(phaseOffsets{1,i}(jumps(j))), ' - ',num2str(phaseOffsets{1,i}(jumps(j+1)-1)),'. '];
            end
        else
            poStr = [num2str(phaseOffsets{1,i}(1)),' & ',num2str(phaseOffsets{1,i}(2)),'. '];
        end
        legendStr{i} = ['Channel for POs: ', poStr];
    end
    legend(lines,legendStr,'Location','best')
    
end

results.paths = paths_ordered;
results.pathlengths = pathlengths_ordered;
results.pathorders = indices;
results.drivPO = poColl;
if length(dvals) > 1
    results.(dv) = dvColl;
end

end

