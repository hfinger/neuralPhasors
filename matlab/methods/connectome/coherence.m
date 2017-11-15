function [ coh ] = coherence( simResult, washOut, windows, fullCoh)
%COHERENCE calculates the coherence of the simulation results
%   simResult            = structure including simulation results,
%                          parameters, Connectivity and Distance matrix
%   washOut              = initial washOut time to remove
%   windows              = time windows between which to calculate mean
%   stimRegionOnly       = if true, only evaluate coherence of each region
%                           with the stimulated region
%
%   Returns:
%       coh              = 3 dimensional coherence matrix

%% calculate coherence

sig = simResult.Y(:,washOut:end);
n = size(sig,1);    
nWindows = size(windows,2);

if ~fullCoh
    
    % calculate coherence between driven nodes and all others
    d = length(simResult.sim.drivPos);
    coh = zeros(d,n,nWindows);
    for i=1:d
        % get phase difference
        diff = bsxfun(@minus, sig(simResult.sim.drivPos(i),:), sig);
        for j=1:nWindows
            % calculate coherence
            coh(i,:,j) = abs(mean(exp(1i*diff(:,windows(1,j):windows(2,j))), 2));
        end
    end

else
    
    % calculate coherence between all pairs of nodes
    coh = zeros(n,n,nWindows);
    for i=1:n
        % get phase difference
        diff = bsxfun(@minus, sig(i,:), sig);
        for j=1:nWindows
            % calculate coherence
            coh(i,:,j) = abs(mean(exp(1i*diff(:,windows(1,j):windows(2,j))), 2));
        end
    end

end

coh = squeeze(mean(coh,3));

end

