function [ SE ] = shanEntrop( sig, nBins, normSE )
%SHANENTROP Function that calculates the shannon entropy of a signal
%
%   Input Parameters:
% 
%       sig - 2-dim array including phase values, with 2. dim = samples/time
%       nBins - number of bins to use to discretize phase values
%       normSE - if true, normalized SE will be returned
%
%   Output:
% 
%       SE - Vector of length N, with N = number of rows of sig

%%

sig = squeeze(sig);
n_areas = size(sig,1);
SE = zeros(n_areas,1);

for i=1:n_areas
    
    distr = histcounts(sig(i,:), linspace(0,2*pi,nBins));
    probs = distr / sum(distr);
    SE(i) = -sum(probs(probs ~= 0) .* log(probs(probs ~= 0)));
    if normSE == 1
        Smax = log(nBins);
        SE(i) = (Smax - SE(i)) / Smax;
    end

end

