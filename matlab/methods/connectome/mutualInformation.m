function [ MI ] = mutualInformation( simResult, washOut, fullMI, nWindows, nBins )
%MUTUALINFORMATION Function that uses the mutInfo function of the
%information theory toolbox to calculate the mutual information between
%rows of the signal stored in simResult. Signal is going to be discretized
%into nBins bins.
%
%   Input Parameters:
%
%       simResult - structure that includes field Y (signal)
%       washOut - initial time to cut off from signal
%       fullMI - if false, MI will only be calculated for driven regions to
%                all other regions
%       nWindows - number of timeWindows over which to calculate MI
%       nBins - number of bins to use to discretize the signal
% 
%   Output:
%
%       MI - 3-dimensional array with the first 2 dimensions being the
%            rows of the signal and the last one being the time windows 

%% extract signal and initialize stuff

sig = simResult.Y(:,washOut:end);
n = size(sig,1);
tmax = size(sig,2);
windowLength = floor(tmax/nWindows);

%% calculate mutual information

% discretize variables
edges = linspace(min(min(sig)), max(max(sig)), nBins);
values = bsxfun(@plus,edges,[edges(2:end),0])/2;
sigDiscrete = discretize(sig,edges,values(1:end-1));

% if not fullMI, calculate mutual information only for driven regions 
if ~fullMI
    
    % get driven regions
    drivPos = simResult.sim.drivPos;
    d = length(drivPos);
    MI = zeros(d,n,nWindows);
    
    % loop over regions
    for i=1:d
        for j=1:n
        % loop over time windows
            for k=1:nWindows
                MI(i,j,k) = mutInfo(sigDiscrete(drivPos(i),(k-1)*windowLength+1:k*windowLength),sigDiscrete(j,(k-1)*windowLength+1:k*windowLength));
            end
        end
    end

else

    MI = zeros(n,n,nWindows);
    % loop over regions
    for i=1:n
        for j=1:n
        % loop over time windows
            for k=1:nWindows
                MI(i,j,k) = mutInfo(sigDiscrete(i,(k-1)*windowLength+1:k*windowLength),sigDiscrete(j,(k-1)*windowLength+1:k*windowLength));
            end
        end
    end

end


end

