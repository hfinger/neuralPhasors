function [ hg ] = plotBrainConnOfSubj( SubjData, ED, driverID, minCoh, plotEntrainmentOnly )
%PLOTBRAINCONNOFSUBJ Plots coherence between different brain regions of
%                    certain subject/run
%   Input Parameters:
%       SubjData - Structure that contains Coherence_full and simResult as
%                  fields
%       ED - Matrix containing euclidean distance between each pair of nodes
%       driverID - Index of the driver to investigate (Scalar)
%       minCoh - Minimum value of coherence for which to plot a connection
%                between two nodes (Scalar)
%       plotEntrainmentOnly - If True, only show coherence between driven
%                             node and all other nodes
%   Output:
%       hg - figure handle

%% Get stimulation strength for each node in the network
drivPos = SubjData.stimPos(driverID);
drivStrength = normpdf(ED(drivPos,:), 0, SubjData.stimRange);
drivStrength = drivStrength - min(drivStrength);
drivConn = (drivStrength/max(drivStrength)) * SubjData.stimScale;

%% Extract coherence data and plot brain connectivity
nDrivers = length(SubjData.stimPos);
Coh = SubjData.Coherence(nDrivers+1:end,nDrivers+1:end);

if plotEntrainmentOnly
    Coh(1:end ~= drivPos,:) = 0;
end

hg = plotBrainConnectivity(drivConn, Coh, minCoh, 1);

end

