function [ SE ] = shannonEntropy( Phase, nBins, tWindows, normSE, targets )
%SHANNOPNENTROPY Function that calculates the shannon entropy for each pair
%of nodes
%
%   Input Parameters:
%       
%       Phase - 2-dim array, with phase signal over time (2.dim)
%       nBins - number of bins to use to discretize the phase signal
%       tWindows - time windows for which to calculate SE
%       normSE - if true, normalized SE will be calculated
%       targets - if not empty, FC will be calculated only for target nodes
% 
%   Output:
%   
%       SE - NxN array with N being the length of the first dimension of
%            simResult.Y or the number of targets

%% get signal

n = size(Phase,1);
m = size(tWindows,2);

%% calculate SE

if isempty(targets)
    targets = 1:n;
end
    
% calculte SE for each pair of nodes
FC = zeros(length(targets),n,m);

for j=1:length(targets)

    % calculate phase difference
    diff = bsxfun(@minus, Phase(targets(j),:), Phase);
    PhaseDiff = abs(mod(diff,2*pi));

    % calculate SE
    for i=1:m
        FC(j,:,i) = shanEntrop(PhaseDiff(:,tWindows(1,i):tWindows(2,i)),nBins,normSE);
    end

end

SE = squeeze(mean(FC,3));

end

