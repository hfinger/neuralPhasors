function [ X ] = addDrivers( X, drivPos, drivRange, drivScale, nodeCoordinates, delay )
%addDrivers adds a number of external drivers to a (connectivity) matrix X
%   X - connectivity or Distance matrix (N x N)
%   drivPos - scalar or vector with positions in X that should be driven
%   drivRange - range of the gaussian to  be used to moddel driver decay
%   drivScale - strength of the drivers
%   ED - (N x N) matrix with the euclidean distances between nodes
%   delay - if True, driver delay of 1 will be added for each node instead
%           of the driver strength (assumes Delay matrix has been passed)
%   Returns new X including the drivers

%% create vector with driver strength for each node generated from a gaussian centered at eeg electrode position nearest to drivPos

% load eeg electrode coordinates and get distances to nodes
electrodeCoordinatesTable = readtable('/net/store/nbp/projects/phasesim/databases/SC_Bastian/surfaces/ca_electrodeLocations/ca03_EEGLocations.txt','Delimiter',' ','ReadVariableNames',false);
electrodeCoordinates = table2array(electrodeCoordinatesTable(:,2:end));
ED = pdist2(electrodeCoordinates, nodeCoordinates);

% find electrode position closest to node to be stimulated and extract
% stimulation strength for each node from gaussian centered around that
% electrode position
nDrivers = size(drivPos',1);
drivConn = zeros(nDrivers,size(X,1) + nDrivers,1);
for i = 1:nDrivers
    if delay
        drivConn(i,nDrivers+1:end) = ones(1,size(X,1));
    else
        [~,pos] = min(ED(:,drivPos(i)));
        drivStrength = normpdf(ED(pos,:), 0, drivRange);
        drivStrength = drivStrength - min(drivStrength);
        drivConn(i,nDrivers+1:end) = (drivStrength/max(drivStrength)) * drivScale;
    end
end

%% include driver in X
X_tmp = zeros(size(X,1) + nDrivers);
X_tmp(1+nDrivers:end,1+nDrivers:end) = X;
X_tmp(:,1:nDrivers) = drivConn';
X = X_tmp;

end

