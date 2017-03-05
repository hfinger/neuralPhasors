function [ SE ] = shanEntrop( coh, nBins, normSE )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

coh = squeeze(coh);
n_areas = size(coh,1);
SE = zeros(n_areas,1);

for i=1:n_areas
    
    distr = histcounts(coh(i,:), linspace(0,2*pi,nBins));
    probs = distr / sum(distr);
    SE(i) = -sum(probs(probs ~= 0) .* log(probs(probs ~= 0)));
    if normSE == 1
        Smax = log(nBins);
        SE(i) = (Smax - SE(i)) / Smax;
    end

end

