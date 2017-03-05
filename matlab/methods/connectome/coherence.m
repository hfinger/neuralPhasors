function [ coh ] = coherence( simResult, washOut, meanCoh, windows, fullCoh)
%COHERENCE calculates the coherence of the simulation results
%   simResult            = structure including simulation results,
%                          parameters, Connectivity and Distance matrix
%   washOut              = initial washOut time to remove
%   meanCoh              = if true, calculate mean Coherence, if false just return pairwise differences 
%   windows              = time windows between which to calculate mean
%   stimRegionOnly       = if true, only evaluate coherence of each region
%                           with the stimulated region
%
%   Returns:
%       coh              = 3 dimensional coherence matrix

sig = simResult.Y(:,washOut:end);
n = size(sig,1);
nWindows = size(windows,2);

if meanCoh
    
    if ~fullCoh
        
        coh = zeros(length(simResult.sim.drivPos),n,nWindows);
        for i=1:length(simResult.sim.drivPos)
            diff = bsxfun(@minus, sig(simResult.sim.drivPos(i),:), sig);
            for j=1:nWindows
                coh(i,:,j) = abs(mean(exp(1i*diff(:,windows(1,j):windows(2,j))), 2));
            end
        end
        
    else
        
    coh = zeros(n,n,nWindows);
    for i=1:n
        diff = bsxfun(@minus, sig(i,:), sig);
        for j=1:nWindows
            coh(i,:,j) = abs(mean(exp(1i*diff(:,windows(1,j):windows(2,j))), 2));
        end
    end
    
    end
    
else
    
    coh = zeros(n,n,size(sig,2));
    
    for i=1:n
        diff = bsxfun(@minus, sig(i,:), sig);
        coh(i,:,:) = mod(diff, 2*pi);
    end
    
end

end

