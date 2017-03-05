function [ coh_mean ] = distToDriver_coh( C, stimPos, maxDist )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

nStims = length(stimPos);
n = min(length(C)-nStims, maxDist);
coh = zeros(nStims,n);

for s = 1:nStims
    idx = zeros(1,size(C,2));
    idx(stimPos) = 1;
    idx(1:nStims) = 1;
    idx = idx == 1;
    target = stimPos(s);
    SWPs = distance_wei(1./C);
    for t = 1:n
        SWPs(:,idx) = inf;
        [ ~, target ] = min(SWPs(target,:));
        if length(target) > 1
            target = target(1);
        end
        coh(s,t) = C(stimPos(s),target);
        idx(target) = 1;
        if sum(idx) == 0
            break
        end
    end
end

coh_mean = mean(coh,1);

end

