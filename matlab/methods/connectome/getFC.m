function [ FC ] = getFC( simResult, FCMeasures, FCWindows, washout )
%GETFC Function that calculates the functional connectivity between nodes
%based on a signal over time
%
%   Input Parameters:
%       
%       simResult - structure that includes field Y, which is a 2-dim array 
%                   with the signal of each node (1.dim) over time (2.dim)
%       FCMeasures - Cell with strings, indicating which measure to use to 
%                    calculate functional connectivity
%       FCWindows - 2-dim array with time windows for which to calculate SE
%                   1. dim = start,end (2 rows), 2.dim = number of windows
%       washout - initial samples to exclude from FC calculation
% 
%   Output:
%
%       FC - NxN array with N = number of rows of sig

%% pre-calculations

% get signal
sig = simResult.Y(:,washout:end);

% get analytic signal
sigHilbert = zeros(size(sig));
for n=1:size(sigHilbert,1)
    sigHilbert(n,:) = hilbert(sig(n,:));
end

% determine target regions for FC calculation
if simResult.sim.fullFC
    targets = 1:size(sigHilbert,1);
else
    targets = simResult.sim.drivPos;
end

% scale time windows
FCWindows = FCWindows / (simResult.sim.dt * simResult.sim.sampling);
    
%% calculate functional connectivity

% loop over all fc measures of interest
FC = cell(1,length(FCMeasures));
for m=1:length(FC)
    
    FCMeasure = FCMeasures{m};
    
    if strcmp(FCMeasure,'PLV_Holger')
        
        % normalize signal with envelope 
        sigPhase = sigHilbert./abs(sigHilbert);
        
        FC_tmp = zeros(length(targets),size(sigHilbert,1),size(FCWindows,2));
        % calculate phase locking value
        for n=1:size(FCWindows,2)
            for k=targets
                FC_tmp(k,:,n) = abs(mean(bsxfun(@times, conj(sigPhase(k,FCWindows(1,n):FCWindows(2,n))), sigPhase(:,FCWindows(1,n):FCWindows(2,n))) , 2));
            end
        end
        FC_tmp = squeeze(mean(FC_tmp,3));
        
    elseif strcmp(FCMeasure,'PLV')
        
        % get phase from analytic signal
        sigPhase = angle(sigHilbert);
        
        FC_tmp = zeros(length(targets),size(sigHilbert,1),size(FCWindows,2));
        % calculate phase locking value
        for n=1:size(FCWindows,2)
            for k=targets
                FC_tmp(k,:,n) = abs(mean(exp(1i*bsxfun(@minus, sigPhase(k,FCWindows(1,n):FCWindows(2,n)), sigPhase(:,FCWindows(1,n):FCWindows(2,n)))) , 2));
            end
        end
        FC_tmp = squeeze(mean(FC_tmp,3));
  
    elseif strcmp(FCMeasure,'phaseDiffHist')
        
        % get phase from analytic signal
        sigPhase = angle(sigHilbert);
        
        FC_tmp = zeros(length(targets),size(sigHilbert,1),size(FCWindows,2),simResult.sim.nBins);
        histEdges = linspace(0, 2*pi, simResult.sim.nBins+1);
        % calculate phase locking value
        for n=1:size(FCWindows,2)
            for k=targets
                phaseDiffs = bsxfun(@minus, sigPhase(k,FCWindows(1,n):FCWindows(2,n)), sigPhase(:,FCWindows(1,n):FCWindows(2,n)));
                phaseDiffs = mod(phaseDiffs, 2*pi);
                tmpPhaseDiffs = num2cell(phaseDiffs,2);
                histResults = cellfun(@(x) histcounts(x,histEdges), tmpPhaseDiffs, 'UniformOutput', false);
                FC_tmp(k,:,n,:) = cell2mat(histResults);
            end
        end
        FC_tmp = squeeze(mean(FC_tmp,3));

    elseif strcmp(FCMeasure,'Coherence') || strcmp(FCMeasure,'Coherency')
        
        FC_tmp = zeros(length(targets),size(sigHilbert,1),size(FCWindows,2));
        for n=1:size(FCWindows,2)
            
            % calculate autospectrum
            autoSpec = mean(conj(sigHilbert(:,FCWindows(1,n):FCWindows(2,n))) .* sigHilbert(:,FCWindows(1,n):FCWindows(2,n)) , 2);

            % calculate crossspectrum
            crossSpec = zeros(length(targets),size(sigHilbert,1));
            for k=targets
                crossSpec(k,:) = mean(bsxfun(@times, conj(sigHilbert(k,FCWindows(1,n):FCWindows(2,n))), sigHilbert(:,FCWindows(1,n):FCWindows(2,n))),2);
            end

            % calculate Coherency
            FC_tmp(:,:,n) = crossSpec ./ sqrt( autoSpec * autoSpec' );

            if strcmp(FCMeasure,'Coherence')
                FC_tmp(:,:,n) = abs(FC_tmp(:,:,n));
            end
            
        end
        FC_tmp = squeeze(mean(FC_tmp,3));
    
    elseif strcmp(FCMeasure,'SE')
        
        % calculate normalized shannon entropy
        phaseSig = angle(sigHilbert);
        FC_tmp = shannonEntropy(phaseSig,simResult.sim.nBins,FCWindows,1,targets);
        
    elseif strcmp(FCMeasure,'MI')
        
        % calculate mutual information
        FC_tmp = mutualInformation(simResult,1,simResult.sim.fullFC,simResult.sim.nWindows,simResult.sim.nBins);
        
    end

    FC{m} = FC_tmp;
    
end

end