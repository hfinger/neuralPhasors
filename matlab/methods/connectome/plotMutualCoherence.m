function [ mutualCoherence, dvSorted ] = plotMutualCoherence( dataStruct, targetPO, dv, pltCoh, drivPosCoh )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(dataStruct);
n = length(fnames);
dataTmp = dataStruct.(fnames{1});
drivPos = dataTmp.drivPos;

mutualCoherence = zeros(n,size(dataTmp.Coherence{1,1}(:,length(drivPos)+1:end),2));
POs = zeros(1,n);
dvColl = zeros(1,n);
for f=1:n
    data = dataStruct.(fnames{f});
    POs(f) = data.drivPO;
    dvColl(f) = data.(dv);
    Coh = data.Coherence{1,1}(:,length(drivPos)+1:end);
    mutualCoherence(f,:) = Coh(drivPos(1),:) .* Coh(drivPos(2),:);
end


uniquePOs = unique(POs);
if strcmp(targetPO,'minmax')
    
    [~,minIdx] = min(drivPosCoh,[],2);
    [~,maxIdx] = max(drivPosCoh,[],2);
    
    uniqueDV = unique(dvColl);
    cohDiff = zeros(length(uniqueDV),size(mutualCoherence,2));
    for d=1:length(uniqueDV)
        minPO = uniquePOs(minIdx(d));
        maxPO = uniquePOs(maxIdx(d));
        idx1 = [POs == minPO] .* [dvColl == uniqueDV(d)];
        idx2 = [POs == maxPO] .* [dvColl == uniqueDV(d)];
        coh1 = mutualCoherence(idx1 == 1,:);
        coh2 = mutualCoherence(idx2 == 1,:);
        cohDiff(d,:) = coh2 - coh1;
    end
    mutualCoherence = cohDiff;
    dvSorted = uniqueDV;
    
else
    
    if sum(uniquePOs == targetPO) == 0
        targetPO = uniquePOs(diff(uniquePOs < targetPO) == -1);
    end
    
    idx1 = POs == targetPO;
    mutualCoherence = mutualCoherence(idx1,:);
    dvColl = dvColl(idx1);
    [dvSorted, idx2] = sort(dvColl);
    mutualCoherence = mutualCoherence(idx2,:);

end


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

