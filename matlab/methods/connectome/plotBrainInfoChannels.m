function [ results ] = plotBrainInfoChannels( paths, dataStruct, k, kPlt, dvs, dvTargets, nodeVals, nodeMin, nodeScale, edgeMin, edgeScale, invCoh )
%PLOTBRAININFOCHANNELS Function that plots a number of paths through the
%connectome specified by paths, with edges based on coherence data in
%dataStruct. Paths will be re-arranged according to coherence data
%(shortest to longest) for each plot.
%   Input Parameters:
%       paths - cell containing sorted path vectors (shortest to longest)
%       dataStruct - structure containing simulation data (including
%                    Coherence)
%       k - number of shortest paths to take into consideration (scalar)
%       kPlt - number of shortest paths to plot in connectome (scalar)
%       dvs - cell with strings indicating for which variables to create
%             plots (max 2 variables)
%       dvTargets - array with one row for each variable in dvs. Specifies
%                   the values of each dv for which to create plots
%                   (if values should only be specified for one variable,
%                   leave the other dvTarget row as zeros)
%       nodeVals - array containing the values used for size and coloring
%                  of the network's nodes. Can be a vector - will then be
%                  used for each plot. Alternatively, use M1xM2xN array,
%                  with M1 being the number of values of dv 1 and M2
%                  the number of values of dv 2
%       nodeMin - Minimum node size (scalar)
%       nodeScale - Scaling of node size with nodeVals (scalar)
%       edgeMin - Minimum value a edge needs to have to be plotted (scalar)
%       edgeScale - Scaling of edge thickness with edgeVals (scalar)
%       invCoh - if true, coherence values will be inverted before using
%                them to find the shortest path between the target nodes
%
%   Output:
%       results - structure, containing the paths, the path
%       activations, the indices of the original paths and the dependent
%       variable values

%% Go through structure and extract relevant data. Then calculate path activations of paths based on Coherence and re-arrange paths from shortest to longest

% initializations
fnames = fieldnames(dataStruct);
n = length(fnames);
paths_ordered = cell(n,k);
pathActivations_ordered = zeros(n,k);
pathInds = zeros(n,k);
dvColl = zeros(n,length(dvs));

% for each structure entry
for f=1:length(fnames)
    
    data = dataStruct.(fnames{f});
    
    % get target nodes
    targets = data.drivPos;
    targets = targets - length(targets);
    
    % get Coherence data
    Coh = data.Coherence{1,1}(length(targets)+1:end,length(targets)+1:end);
    if invCoh
        Coh = 1./ Coh;
    end
    
    % get dependent variables
    for d=1:length(dvs)
        dvColl(f,d) = data.(dvs{d});
    end
    
    % compute path activations of shortest paths based on Coherence matrix
    pathActivations = zeros(1,k);
    for i=1:k
        for j=1:length(paths{1,i})-1
            pathActivations(i) = pathActivations(i) + Coh(paths{1,i}(j),paths{1,i}(j+1));
        end
        if ~invCoh
            pathActivations(i) = pathActivations(i) / (length(paths{1,i})-1);
        end
    end

    % order the path activations and get path indices
    [pathActivations_ordered(f,:),idx] = sort(pathActivations);
    if ~invCoh
        idx = fliplr(idx);
        pathActivations_ordered(f,:) = pathActivations_ordered(f,idx);
    end
    
    % order paths according to path indices
    for i=1:length(idx)
        paths_ordered{f,i} = paths{1,idx(i)};
    end
    pathInds(f,:) = idx;

end

if invCoh
    pathActivations_ordered = 1./pathActivations_ordered;
end

%% Go through all combinations of dependent variable values and plot paths in connectome

% get relevant values of first dependent variable
dv1Vals = unique(dvColl(:,1));
if ~isempty(dvTargets) && sum(dvTargets(1,:)) ~= 0
    dv1Vals_tmp = dvTargets(1,:);
    for i=1:length(dv1Vals_tmp)
        if sum(dv1Vals == dv1Vals_tmp(i)) == 0
            dv1Vals_tmp(i) = dv1Vals(diff(dv1Vals < dv1Vals_tmp(i)) == -1);
        end
    end
    dv1Vals = unique(dv1Vals_tmp);
end

% set plotting parameters
nodeRange = [min(min(min(nodeVals))),max(max(max(nodeVals)))];
surfaceVisibility = 0.1;
edgeAlphas = linspace(0.2,1.0,kPlt);

if isempty(nodeVals)
    nodeVals = zeros(1,size(Coh,1));
    nodeVals(targets) = 1;
end
nodeScale = nodeScale/max(max(max(nodeVals)));
nodeMinSize = nodeMin - min(min(min(nodeVals)))*nodeScale;

% loop over first dependent variable
for d1=1:length(dv1Vals)
    
    % get data where dv1 has specified value
    idx1 = dvColl(:,1) == dv1Vals(d1);
    dvColl_tmp = dvColl(idx1,:);
    paths_ordered_tmp = paths_ordered(idx1,:);
    pathInds_tmp = pathInds(idx1,:);
    fnames_tmp = fnames(idx1);
    
    % get relevant values of second dependent variable
    dv2Vals = unique(dvColl(:,2));
    if ~isempty(dvTargets) && sum(dvTargets(2,:)) ~= 0
        dv2Vals_tmp = dvTargets(2,:);
        for i=1:length(dv2Vals_tmp)
            if sum(dv2Vals == dv2Vals_tmp(i)) == 0
                dv2Vals_tmp(i) = dv2Vals(diff(dv2Vals < dv2Vals_tmp(i)) == -1);
            end
        end
        dv2Vals = unique(dv2Vals_tmp);
    end
    
    if length(size(nodeVals)) < 3
       nodeVals = reshape(nodeVals,1,1,length(nodeVals));
       nodeVals = repmat(nodeVals,length(dv1Vals),length(dv2Vals),1);
    end
    
    % loop over second dependent variable
    for d2=1:length(dv2Vals)
        
        % get data where dv2 has specified value
        idx2 = dvColl_tmp(:,2) == dv2Vals(d2);
        paths_ordered_tmp2 = paths_ordered_tmp(idx2,:);
        pathInds_tmp2 = pathInds_tmp(idx2,:);
        fnames_tmp2 = fnames_tmp(idx2);
        
        % loop over number of paths to plot
        edgeVals = zeros(size(Coh,1),size(Coh,2),kPlt);
        for path=kPlt:-1:1

            % get coherence of simulation to plot
            Coh = dataStruct.(fnames_tmp2{1}).Coherence{1,1}(length(targets)+1:end,length(targets)+1:end);
            % get edge values for path plotting
            for j=1:length(paths_ordered_tmp2{1,pathInds_tmp2(path)})-1;
                i1 = paths_ordered_tmp2{1,pathInds_tmp2(path)}(j);
                i2 = paths_ordered_tmp2{1,pathInds_tmp2(path)}(j+1);
                %coh = Coh(i1,i2);  
                edgeVals(i1,i2,kPlt-path+1) = (kPlt-path+1)*edgeScale;
            end
        end
        
        % plot kPlt paths in connectome
        figure()
        plotBrainConnectivity(squeeze(nodeVals(d1,d2,:))',edgeVals,nodeMinSize,nodeScale,nodeRange,'r',edgeMin,edgeAlphas,surfaceVisibility,0,1);
        if kPlt == 1
            set(gcf,'name',['Shortest path for ',dvs{1},' = ',num2str(dv1Vals(d1)),' & ',dvs{2},' = ',num2str(dv2Vals(d2))])
        else
            set(gcf,'name',[num2str(kPlt),'Shortest paths for ',dvs{1},' = ',num2str(dv1Vals(d1)),' & ',dvs{2},' = ',num2str(dv2Vals(d2))])    
        end
        
    end
    
end

%% create results structure

results.paths = paths_ordered;
results.pathActivations = pathActivations_ordered;
results.pathorders = pathInds;
results.(dvs{1}) = dvColl(:,1);
results.(dvs{2}) = dvColl(:,2);

end

