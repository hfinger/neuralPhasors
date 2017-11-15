function [ bestTrial, bestTrialInd, dvColl ] = getBestTrial( struct, dv )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(struct);
n = length(fnames);
dvColl = zeros(n,1);

for f=1:n
   dvColl(f) = struct.(fnames{f}).(dv);
end
[~,ind] = max(dvColl,[],1);

bestTrialInd = fnames{ind(1)};
bestTrial = struct.(fnames{ind(1)});

end

