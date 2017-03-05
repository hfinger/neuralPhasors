function [ bestTrial, bestTrialInd, dvColl ] = getBestTrial( struct, dv )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fnames = fieldnames(struct);
n = length(fnames);
m = length(struct.(fnames{1}).(dv));
dvColl = zeros(n,m);

for f=1:n
   data = struct.(fnames{f}).(dv);
   for i=1:m
        dvColl(f,i) = data{1,i};
   end
end
[~,ind] = max(dvColl,[],1);

for j=1:length(ind)
    bestTrialInd{j} = fnames{ind(j)};
    bestTrial{j} = struct.(fnames{ind(j)});
end

end

