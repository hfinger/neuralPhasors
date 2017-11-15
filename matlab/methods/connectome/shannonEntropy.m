function [ SE ] = shannonEntropy( simResult, nBins, tWindows, normSE, fullFC )
%SHANNOPNENTROPY Function that calculates the shannon entropy for each pair
%of nodes
%
%   Input Parameters:
%       
%       simResult - structure, including field Y with signal over time
%       nBins - number of bins to use to discretize the phase signal
%       tWindows - time windows for which to calculate SE
%       normSE - if true, normalized SE will be calculated
% 
%   Output:
%   
%       SE - NxN array with N being the length of the first dimension of
%            simResult.Y

%% get signal

Phase = simResult.Y;
n = size(Phase,1);
m = size(tWindows,2);

%% calculate SE

if fullFC
    
    % calculte SE for each pair of nodes
    clear simResult
    FC = zeros(n,n,m);
    
    for j=1:n
        
        % calculate phase difference
        diff = bsxfun(@minus, Phase(j,:), Phase);
        PhaseDiff = mod(diff, 2*pi);
        
        % calculate SE
        for i=1:m
            FC(j,:,i) = shanEntrop(PhaseDiff(:,tWindows(1,i):tWindows(2,i)),nBins,normSE);
        end
        
    end
    
else
    
    % calculate SE for each driven node
    drivPos = simResult.sim.drivPos;
    clear simResult
    FC = zeros(length(drivPos),n,m);
    
    for j=1:length(drivPos)
        
        % calculate phase difference
        diff = bsxfun(@minus, Phase(drivPos(j),:), Phase);
        PhaseDiff = mod(diff, 2*pi);
        
        %c alculate SE
        for i=1:m
            FC(j,:,i) = shanEntrop(PhaseDiff(:,tWindows(1,i):tWindows(2,i)),nBins,normSE);
        end
        
    end
    
end

SE = squeeze(mean(FC,3));

end

