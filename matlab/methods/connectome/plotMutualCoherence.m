function [ mutualCoherence, dvSorted ] = plotMutualCoherence( dataStruct, targetPO, dv, pltCoh )
%PLOTMUTUALCOHERENCE Function that plots the coherence of all network nodes
% to two driven regions for a certain phase offset over all unique values
% of a given dependent variable
% 
%   Input Parameters:
%
%       dataStruct - structure that contains 1 field for each simulation
%       targetPO - Either the phase offset value, for which to evaluate the mutual
%                  coherence (scalar), or 'minmax' (string). In the second
%                  case, the difference in mutual coherence between the
%                  PO for which the coherence between the driven regions is
%                  maximal is evaluated.
%       dv - name of the dependent variable (string)
%       pltCoh - if false, no plot is created
%
%   Output:
%
%       mutualCoherence - N by M array of the mutual coherence, where N is
%       the number of unique values of dv and M is the number of network
%       nodes
%       dvSorted - vector of length N with the sorted unique values of dv

%% extract coherence, phase offset and dv from struct

fnames = fieldnames(dataStruct);
n = length(fnames);
dataTmp = dataStruct.(fnames{1});
drivPos = dataTmp.drivPos;

mutualCoherence = zeros(n,size(dataTmp.Coherence{1,1}(:,length(drivPos)+1:end),2));
POs = zeros(1,n);
dvColl = zeros(1,n);
drivPosCoh = zeros(1,n);
for f=1:n
    data = dataStruct.(fnames{f});
    POs(f) = data.drivPO;
    dvColl(f) = data.(dv);
    drivPosCoh(f) = data.Coherence{1,1}(drivPos(1),drivPos(2));
    Coh = data.Coherence{1,1}(:,length(drivPos)+1:end);
    mutualCoherence(f,:) = Coh(drivPos(1),:) .* Coh(drivPos(2),:);
end

%% calculate mutual coherence for target phase offset and all values of dv

% if difference between minimum and maximum coherence phase offset should
% be evaluated
if strcmp(targetPO,'minmax')
    
    % get unique values of dv and loop over them
    uniqueDV = unique(dvColl);
    cohDiff = zeros(length(uniqueDV),size(mutualCoherence,2));
    for d=1:length(uniqueDV)
        
        % extract fields where dv takes target value d
        idx = dvColl == uniqueDV(d);
        drivPosCoh_tmp = drivPosCoh(idx);
        mutualCoherence_tmp = mutualCoherence(idx,:);
        
        % get index of phase offset at min and max coherence
        [~,minIdx] = min(drivPosCoh_tmp);
        [~,maxIdx] = max(drivPosCoh_tmp);
        
        % get mutual coherence difference
        cohMin = mutualCoherence_tmp(minIdx,:);
        cohMax = mutualCoherence_tmp(maxIdx,:);
        cohDiff(d,:) = cohMax - cohMin;
        
    end
    mutualCoherence = cohDiff;
    dvSorted = uniqueDV;

% else get mutual coherence at target phase offset
else
    
    % get unique phase offsets
    uniquePOs = unique(POs);

    % look for unique PO value nearest to target PO
    if sum(uniquePOs == targetPO) == 0
        targetPO = uniquePOs(diff(uniquePOs < targetPO) == -1);
    end
    
    % get mutual coherence at target PO and sort it after dv
    idx1 = POs == targetPO;
    mutualCoherence_tmp = mutualCoherence(idx1,:);
    dvColl_tmp = dvColl(idx1);
    [dvSorted, idx2] = sort(dvColl_tmp);
    mutualCoherence = mutualCoherence_tmp(idx2,:);

end

%% Plotting

if pltCoh
    figure()
    imagesc(mutualCoherence)
    yticklabels = cell(1,length(dvSorted));
    for i=1:length(dvSorted)
        yticklabels{i} = num2str(dvSorted(i));
    end
    set(gca, 'YTick', 1:length(dvSorted),'YTickLabel',yticklabels)
    xlabel('Network Nodes')
    ylabel(dv)
    title(['Coherence with driven nodes for phase offset = ', num2str(targetPO)])
end

end

