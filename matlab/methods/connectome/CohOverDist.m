function [ CohDist, Targets ] = CohOverDist( Coh, SC, stimPos, stimNodes, maxDist )
%COHOVERDIST Function that evaluates coherence over distance as measured by
%            shortest weighted paths based on the structural connectivity
%            matrix
%   Input Parameters:
%       Coh - N x N functional connectivity matrix with N = number of Nodes
%       SC - N x N structural connectivity matrix
%       stimPos - stimulated node for which to evaluate coherence over
%                 distance
%       stimNodes - all nodes stimulated during simulation (vector)
%       maxDist - maximum distance for which to evaluate coherence (scalar)
%   
%   Output:
%       CohDist - Vector of length maxDist containing coherence values for
%                 each distance
%       Targets - Vector of length maxDist containing node indices for
%                 which coherence was evaluated

%%

% get shortest weighted paths for each pair of nodes
SWPs = distance_wei(1./SC);

% which nodes to include in evaluation (1 means to exclude node)
ind = zeros(1,length(SC));
ind(stimNodes) = 1;

CohDist = zeros(1,maxDist);
Targets = uint8(zeros(1,maxDist));

% extract coherence for each step away from stimulated node
for d=1:maxDist
    
    SWPs(stimPos,ind == 1) = inf;
    [~,target] = min(SWPs(stimPos,:));
    CohDist(d) = Coh(stimPos,target);
    Targets(d) = target;
    
    ind(target) = 1;
    if sum(ind == 0) == 0
        break
    end
    
end

end

